"""Unified genotype likelihood: read depth + junction support + BAF.

For a candidate DEL/DUP with genotypes g in {0/0, 0/1, 1/1} (diploid contig)
or {0, 1} (haploid), the posterior combines:

- RD:   Gaussian likelihood of the per-bin copy ratio around the expected
        ratio of g, variance scaling with the mean (as in the HMM), total
        depth weight capped so megabase events don't saturate;
- PE/SR: Poisson likelihood of the junction-supporting read count with rate
        proportional to the (haploid) dosage of the variant allele;
- BAF:  beta(-mixture) likelihood of clear-het VAFs inside the event under
        the copy state implied by g.

This is deliberately a *generative, interpretable* model: each term is also
exported as a feature to the downstream calibrated filter.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln
from scipy.stats import beta as beta_dist

from ..config import GenotypeConfig
from ..variants import CNVCall

logger = logging.getLogger(__name__)

_MIN_RATIO_VAR = 0.05
_MAX_DEPTH_BINS = 200      # cap on effective depth-informative bins
_MAX_BAF_SITES = 50        # cap on effective BAF sites
_JUNCTION_EFFICIENCY = 0.12  # detectable junction reads per unit of haploid coverage


def expected_ratio(alt_dosage: int, svtype: str, ploidy: int) -> float:
    """Expected copy ratio (diploid-normalized) given alt allele dosage."""
    cn = ploidy + alt_dosage * (1 if svtype == "DUP" else -1)
    return max(cn, 0) / 2.0


def _ratio_loglik(ratio: np.ndarray, expected: float, sigma0: float) -> float:
    var = sigma0**2 * max(expected, _MIN_RATIO_VAR)
    ll = -0.5 * np.log(2 * np.pi * var) - 0.5 * (ratio - expected) ** 2 / var
    if len(ll) > _MAX_DEPTH_BINS:
        ll = ll * (_MAX_DEPTH_BINS / len(ll))
    return float(np.sum(ll))


def _poisson_loglik(count: int, rate: float) -> float:
    rate = max(rate, 1e-6)
    return float(count * np.log(rate) - rate - gammaln(count + 1))


def _baf_loglik(
    vaf: np.ndarray, alt_dosage: int, svtype: str, ploidy: int, mismap: float = 0.01,
    concentration: float = 30.0,
) -> float:
    """Log-likelihood of clear-het VAFs (0.15-0.85) under the implied CN."""
    if len(vaf) == 0:
        return 0.0
    weight = min(1.0, _MAX_BAF_SITES / len(vaf))
    cn = ploidy + alt_dosage * (1 if svtype == "DUP" else -1)
    if cn <= 1:
        return float(np.log(mismap) * len(vaf) * weight)
    means = [i / cn for i in range(1, cn)] or [0.5]
    v = np.clip(vaf, 1e-3, 1 - 1e-3)
    parts = np.stack(
        [beta_dist.logpdf(v, m * concentration, (1 - m) * concentration) for m in means]
    )
    m = parts.max(axis=0)
    site_ll = m + np.log(np.exp(parts - m).mean(axis=0))
    return float(site_ll.sum() * weight)


@dataclass
class ContigData:
    """Per-contig evidence arrays the genotyper slices from."""

    ratio: np.ndarray        # corrected copy ratio per bin
    mask: np.ndarray         # excluded bins
    baf_pos: np.ndarray
    baf_vaf: np.ndarray
    bin_size: int
    ploidy: int
    mean_depth: float        # diploid-normalized units of the ratio denominator
    sigma0: float


class Genotyper:
    def __init__(self, config: GenotypeConfig, max_cn: int = 5):
        self.config = config
        self.max_cn = max_cn

    def _genotypes(self, ploidy: int) -> list[tuple[tuple[int, int], int]]:
        """(GT tuple, alt dosage) per genotype, with priors aligned."""
        if ploidy == 1:
            return [((0,) * 2, 0), ((1, 1), 1)]  # report haploid as 0/0 / 1/1
        return [((0, 0), 0), ((0, 1), 1), ((1, 1), 2)]

    def _priors(self, ploidy: int) -> np.ndarray:
        c = self.config
        if ploidy == 1:
            p = np.array([1 - c.het_prior, c.het_prior])
        else:
            p = np.array([1 - c.het_prior - c.hom_prior, c.het_prior, c.hom_prior])
        return np.log(p)

    def genotype(self, call: CNVCall, data: ContigData) -> None:
        """Assign GT/GQ and integer CN/CNQ to the call (in place)."""
        b0 = call.start // data.bin_size
        b1 = max(-(-call.end // data.bin_size), b0 + 1)
        sel_mask = data.mask[b0:b1]
        ratio = data.ratio[b0:b1][~sel_mask]

        sel_baf = (data.baf_pos >= call.start) & (data.baf_pos < call.end)
        vaf_all = data.baf_vaf[sel_baf]
        hets = vaf_all[(vaf_all > 0.15) & (vaf_all < 0.85)]
        call.baf_n = int(len(hets))
        if len(hets):
            call.baf_dev = float(np.mean(np.abs(hets - 0.5)))

        support = call.pe_support + call.sr_support
        haploid_cov = data.mean_depth / 2.0

        genotypes = self._genotypes(data.ploidy)
        logp = self._priors(data.ploidy).copy()
        for i, (_gt, dosage) in enumerate(genotypes):
            if len(ratio):
                logp[i] += _ratio_loglik(
                    ratio, expected_ratio(dosage, call.svtype, data.ploidy), data.sigma0
                )
            rate = (
                self.config.mismap_rate * data.mean_depth
                if dosage == 0
                else _JUNCTION_EFFICIENCY * haploid_cov * dosage
            )
            # junction term only informs when junction evidence was collected
            if call.source != "DEPTH" or support > 0:
                logp[i] += _poisson_loglik(support, rate)
            logp[i] += self._baf_loglik_safe(hets, dosage, call.svtype, data.ploidy)

        logp -= logp.max()
        post = np.exp(logp)
        post /= post.sum()
        best = int(np.argmax(post))
        call.gt = genotypes[best][0]
        call.gq = int(min(99, -10 * np.log10(max(1 - post[best], 1e-10))))

        # integer CN from depth alone on the final coordinates
        if len(ratio):
            states = np.arange(self.max_cn + 1)
            var = data.sigma0**2 * np.maximum(states / 2.0, _MIN_RATIO_VAR)
            ll = (
                -0.5 * np.log(2 * np.pi * var)[None, :]
                - 0.5 * (ratio[:, None] - (states / 2.0)[None, :]) ** 2 / var[None, :]
            )
            if len(ratio) > _MAX_DEPTH_BINS:
                ll = ll * (_MAX_DEPTH_BINS / len(ratio))
            tot = ll.sum(axis=0)
            tot -= tot.max()
            cn_post = np.exp(tot)
            cn_post /= cn_post.sum()
            call.cn = int(np.argmax(cn_post))
            call.cn_posterior = float(cn_post[call.cn])
            call.mean_ratio = float(ratio.mean())
            call.n_bins = int(b1 - b0)
            call.masked_frac = float(sel_mask.mean())

    def _baf_loglik_safe(
        self, hets: np.ndarray, dosage: int, svtype: str, ploidy: int
    ) -> float:
        try:
            return _baf_loglik(hets, dosage, svtype, ploidy)
        except Exception:  # pragma: no cover - defensive: BAF must never kill a call
            logger.exception("BAF likelihood failed; ignoring BAF term")
            return 0.0
