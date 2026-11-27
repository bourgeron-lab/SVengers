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
from scipy.stats import poisson

from ..config import GenotypeConfig
from ..variants import CNVCall

logger = logging.getLogger(__name__)

_MIN_RATIO_VAR = 0.05
_MAX_DEPTH_BINS = 200      # cap on effective depth-informative bins
_MAX_BAF_SITES = 50        # cap on effective BAF sites
_FLANK_BINS = 5            # bins each side giving an edge bin's copy ratio outside the event
_JUNCTION_EFFICIENCY = 0.12  # detectable junction reads per unit of haploid coverage
_HET_OVERDISPERSION = 2.0   # local het density varies ~2x around the contig rate


def expected_ratio(alt_dosage: int, svtype: str, ploidy: int) -> float:
    """Expected copy ratio (diploid-normalized) given alt allele dosage."""
    cn = ploidy + alt_dosage * (1 if svtype == "DUP" else -1)
    return max(cn, 0) / 2.0


def _ratio_loglik(ratio: np.ndarray, expected: float, sigma0: float,
                  extra_sd: np.ndarray | None = None) -> float:
    var = sigma0**2 * max(expected, _MIN_RATIO_VAR)
    if extra_sd is not None:
        var = var + (extra_sd * expected) ** 2
    ll = -0.5 * np.log(2 * np.pi * var) - 0.5 * (ratio - expected) ** 2 / var
    if len(ll) > _MAX_DEPTH_BINS:
        ll = ll * (_MAX_DEPTH_BINS / len(ll))
    return float(np.sum(ll))


def _poisson_loglik(count: int, rate: float) -> float:
    rate = max(rate, 1e-6)
    return float(count * np.log(rate) - rate - gammaln(count + 1))


def _nb_loglik(count: int, mean: float, alpha: float) -> float:
    """Negative-binomial log-likelihood, variance = mean + alpha*mean^2.

    alpha -> 0 recovers Poisson; alpha > 0 flattens the het-vs-hom contrast a
    raw read count can support, which is what junction-read counts warrant.
    """
    if alpha <= 0:
        return _poisson_loglik(count, mean)
    mean = max(mean, 1e-6)
    r = 1.0 / alpha
    p = r / (r + mean)
    return float(
        gammaln(count + r) - gammaln(r) - gammaln(count + 1)
        + r * np.log(p) + count * np.log1p(-p)
    )


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


def _flank_ratio(data: "ContigData", b0: int, b1: int) -> float:
    """Median copy ratio of the unmasked bins beside [b0, b1) (the contig's
    baseline when none)."""
    lo, hi = max(b0 - _FLANK_BINS, 0), min(b1 + _FLANK_BINS, len(data.ratio))
    side = np.r_[lo:b0, b1:hi]
    side = side[~data.mask[side]]
    return float(np.median(data.ratio[side])) if len(side) else data.ploidy / 2.0


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
    het_rate: float | None = None  # clear hets per bp on the contig; None: no paralog check
    extra_sd: np.ndarray | None = None  # per-bin relative population spread (PoN variance mode)


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
        xsd = data.extra_sd[b0:b1][~sel_mask] if data.extra_sd is not None else None

        sel_baf = (data.baf_pos >= call.start) & (data.baf_pos < call.end)
        vaf_all = data.baf_vaf[sel_baf]
        hets = vaf_all[(vaf_all > 0.15) & (vaf_all < 0.85)]
        call.baf_n = int(len(hets))
        if len(hets):
            call.baf_dev = float(np.mean(np.abs(hets - 0.5)))
        if self._paralog_dense(len(hets), call.end - call.start, data.het_rate):
            hets = hets[:0]

        support = call.pe_support + call.sr_support + call.cig_support
        haploid_cov = data.mean_depth / 2.0

        genotypes = self._genotypes(data.ploidy)
        logp = self._priors(data.ploidy).copy()
        # sub-bin dilution: a 100 bp event only shifts its 1 kb bin's depth by
        # ~10%, so the expected ratio must be pulled toward 1 accordingly —
        # otherwise the depth term vetoes every junction-supported small event
        span = (b1 - b0) * data.bin_size
        dilution = min(1.0, call.size / span) if span else 1.0
        for i, (_gt, dosage) in enumerate(genotypes):
            if len(ratio):
                expected = expected_ratio(dosage, call.svtype, data.ploidy)
                baseline = data.ploidy / 2.0
                diluted = baseline + (expected - baseline) * dilution
                logp[i] += _ratio_loglik(ratio, diluted, data.sigma0, xsd)
            rate = (
                self.config.mismap_rate * data.mean_depth
                if dosage == 0
                else _JUNCTION_EFFICIENCY * haploid_cov * dosage
            )
            # junction term only informs when junction evidence was collected
            if call.source != "DEPTH" or support > 0:
                logp[i] += _nb_loglik(
                    support, rate, self.config.junction_overdispersion
                )
            logp[i] += self._baf_loglik_safe(hets, dosage, call.svtype, data.ploidy)

        # how firmly the data alone (flat prior) rule out the event: a parent
        # with no junction read at a small deletion is 0/0 GQ ~45 by the
        # het prior only, which made inherited polymorphisms missed in both
        # parents look like confirmed de novo
        lik = logp - self._priors(data.ploidy)
        lik -= lik.max()
        p_flat = np.exp(lik)
        p_flat /= p_flat.sum()
        call.rq = int(min(99, -10 * np.log10(max(1 - p_flat[0], 1e-10))))

        logp -= logp.max()
        post = np.exp(logp)
        post /= post.sum()
        best = int(np.argmax(post))
        call.gt = genotypes[best][0]
        call.gq = int(min(99, -10 * np.log10(max(1 - post[best], 1e-10))))
        # conditional on the call being real: the genotype prior (het 1e-4)
        # outweighs the few junction reads of a small event on a 1 kb-depth
        # genome, so a call the filter passes can come out 0/0 here (~1,200
        # per genome); its GT must then be the best non-reference genotype
        nr = logp[1:] - logp[1:].max()
        p_nr = np.exp(nr)
        p_nr /= p_nr.sum()
        k = int(np.argmax(p_nr))
        call.gt_alt = genotypes[1 + k][0]
        call.gq_alt = int(min(99, -10 * np.log10(max(1 - p_nr[k], 1e-10))))

        # integer CN from depth alone on the final coordinates. An edge bin
        # holds the flanks' copy number outside the event: a 1.7 kb CN0 inside
        # a het deletion, read on whole bins, came out CN1 (and a het
        # deletion that size CN2). Below a bin the depth cannot tell the
        # states apart, and such events keep whole-bin reading
        if len(ratio):
            states = np.arange(self.max_cn + 1)
            if call.size >= data.bin_size:
                idx = np.arange(b0, b1)[~sel_mask]
                cover = (np.minimum((idx + 1) * data.bin_size, call.end)
                         - np.maximum(idx * data.bin_size, call.start))
                frac = np.clip(cover / data.bin_size, 0.0, 1.0)[:, None]
                base = _flank_ratio(data, b0, b1)
            else:
                frac, base = np.ones((len(ratio), 1)), 0.0
            expected = base * (1 - frac) + (states / 2.0)[None, :] * frac
            var = data.sigma0**2 * np.maximum(expected, _MIN_RATIO_VAR)
            if xsd is not None:
                var = var + (xsd[:, None] * expected) ** 2
            ll = -0.5 * np.log(2 * np.pi * var) - 0.5 * (ratio[:, None] - expected) ** 2 / var
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

    def _paralog_dense(self, n_het: int, length: int, het_rate: float | None) -> bool:
        """More hets than a diploid stretch of this length yields: they come
        from reads of a paralog mapped here, present whatever the copy number,
        so they say nothing about the genotype (a common deletion's carriers
        would otherwise each pay log(mismap) per het and fall to 0/0)."""
        p = self.config.baf_paralog_p
        if not het_rate or p <= 0 or n_het == 0:
            return False
        expected = _HET_OVERDISPERSION * het_rate * max(length, 1)
        return float(poisson.sf(n_het - 1, expected)) < p

    def _baf_loglik_safe(
        self, hets: np.ndarray, dosage: int, svtype: str, ploidy: int
    ) -> float:
        try:
            return _baf_loglik(hets, dosage, svtype, ploidy)
        except Exception:  # pragma: no cover - defensive: BAF must never kill a call
            logger.exception("BAF likelihood failed; ignoring BAF term")
            return 0.0
