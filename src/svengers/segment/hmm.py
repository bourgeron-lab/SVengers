"""CN-state HMM segmentation of the corrected depth track.

States are integer copy numbers 0..max_cn on a contig of baseline ploidy 1 or 2.
Emissions combine two channels:

- depth: copy ratio r_i = corrected_depth_i / diploid_mean, Gaussian with
  mean c/2 and variance scaling with the mean (Poisson-like), sigma0 fit from
  the data (v1; a negative-binomial count model on 100 bp raw counts is the
  planned upgrade once the self-contained scan is the default path);
- BAF: per-bin log-likelihood of clear-het VAFs under each CN (optional).

Masked bins are uninformative (equal emission under every state), so the chain
persists through gaps instead of breaking segments.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

try:
    from numba import njit
except ImportError:  # pragma: no cover - numba is a hard dep, fallback for docs builds
    def njit(*args, **kwargs):
        def wrap(f):
            return f
        return wrap if not (args and callable(args[0])) else args[0]

_MIN_RATIO_VAR = 0.05  # variance floor multiplier (as a fraction of expected ratio)


@dataclass(frozen=True)
class Segment:
    start_bin: int
    end_bin: int          # exclusive
    cn: int
    posterior: float      # mean per-bin posterior of the called state
    mean_ratio: float     # mean copy ratio over unmasked bins
    n_bins: int
    masked_frac: float
    trim_left: int = 0    # masked bins trimmed off the start (breakpoint uncertainty)
    trim_right: int = 0   # masked bins trimmed off the end


# ------------------------------------------------------------------ core HMM (numba)


@njit(cache=True)
def _viterbi(log_emit, log_stay, log_switch, log_pi):  # pragma: no cover - numba
    n, k = log_emit.shape
    delta = np.empty((n, k))
    psi = np.zeros((n, k), dtype=np.int32)
    delta[0] = log_pi + log_emit[0]
    for i in range(1, n):
        for s in range(k):
            best, arg = -1e300, 0
            for p in range(k):
                v = delta[i - 1, p] + (log_stay if p == s else log_switch)
                if v > best:
                    best, arg = v, p
            delta[i, s] = best + log_emit[i, s]
            psi[i, s] = arg
    path = np.empty(n, dtype=np.int32)
    path[-1] = np.argmax(delta[-1])
    for i in range(n - 2, -1, -1):
        path[i] = psi[i + 1, path[i + 1]]
    return path


@njit(cache=True)
def _forward_backward(log_emit, log_stay, log_switch, log_pi):  # pragma: no cover
    n, k = log_emit.shape
    log_alpha = np.empty((n, k))
    log_beta = np.zeros((n, k))
    log_alpha[0] = log_pi + log_emit[0]
    for i in range(1, n):
        for s in range(k):
            m = -1e300
            for p in range(k):
                v = log_alpha[i - 1, p] + (log_stay if p == s else log_switch)
                if v > m:
                    m = v
            acc = 0.0
            for p in range(k):
                acc += np.exp(log_alpha[i - 1, p] + (log_stay if p == s else log_switch) - m)
            log_alpha[i, s] = m + np.log(acc) + log_emit[i, s]
    for i in range(n - 2, -1, -1):
        for s in range(k):
            m = -1e300
            for p in range(k):
                v = (log_stay if p == s else log_switch) + log_emit[i + 1, p] + log_beta[i + 1, p]
                if v > m:
                    m = v
            acc = 0.0
            for p in range(k):
                acc += np.exp(
                    (log_stay if p == s else log_switch)
                    + log_emit[i + 1, p]
                    + log_beta[i + 1, p]
                    - m
                )
            log_beta[i, s] = m + np.log(acc)
    gamma = log_alpha + log_beta
    for i in range(n):
        m = gamma[i].max()
        gamma[i] = np.exp(gamma[i] - m)
        gamma[i] /= gamma[i].sum()
    return gamma


# -------------------------------------------------------------------------- emissions


def depth_log_emissions(
    ratio: np.ndarray,
    mask: np.ndarray,
    ploidy: int,
    max_cn: int,
    sigma0: float,
    extra_sd: np.ndarray | None = None,
) -> np.ndarray:
    """(n_bins, n_states) Gaussian log-likelihood of the copy ratio.

    ``extra_sd``: per-bin population spread of the ratio relative to its
    centre (PoN variance mode), added as multiplicative noise: variance
    sigma0^2 * max(e, floor) + (extra_sd * e)^2 at expected ratio e.

    Expected ratio for state c is c/2 regardless of contig ploidy (the ratio is
    normalized to the DIPLOID autosomal mean; a hemizygous X has baseline
    CN1 -> ratio 0.5). ``sigma0`` is the ratio sd at CN2.
    """
    states = np.arange(max_cn + 1)
    expected = states / 2.0
    # variance grows with the mean (Poisson-like), floored for CN0
    var = sigma0**2 * np.maximum(expected, _MIN_RATIO_VAR)
    if extra_sd is not None:
        var = var[None, :] + (extra_sd[:, None] * expected[None, :]) ** 2
    else:
        var = var[None, :]
    log_norm = -0.5 * np.log(2 * np.pi * var)
    emit = log_norm - 0.5 * (ratio[:, None] - expected[None, :]) ** 2 / var
    emit[mask] = 0.0  # uninformative bins
    return emit


def fit_sigma0(ratio: np.ndarray, mask: np.ndarray, ploidy: int = 2) -> float:
    """Robust sd of the copy ratio, expressed at CN2 (the emission model's unit).

    Fitted around the contig's own baseline (ploidy / 2) and scaled by
    1/sqrt(baseline), the Poisson-like variance of the emissions: a male chrX
    sits at 0.5, and a 0.7-1.3 window kept only its noisy tail (~170 of 147k
    bins, sigma0 0.02 instead of 0.08), so every dip of a single bin became a
    hemizygous deletion.
    """
    if ploidy <= 0:
        return 0.15
    base = ploidy / 2.0
    vals = ratio[~mask]
    vals = vals[(vals > 0.7 * base) & (vals < 1.3 * base)]
    if len(vals) < 100:
        return 0.15
    mad = float(np.median(np.abs(vals - np.median(vals))) * 1.4826)
    return mad / np.sqrt(base) or 0.15


def baf_log_emissions(
    baf_pos: np.ndarray,
    baf_vaf: np.ndarray,
    n_bins: int,
    bin_size: int,
    max_cn: int,
    concentration: float = 30.0,
    het_lo: float = 0.15,
    het_hi: float = 0.85,
    mismap: float = 0.01,
    max_sites_weight: float = 3.0,
) -> np.ndarray:
    """(n_bins, n_states) log-likelihood of observed clear-het VAFs under each CN.

    Expected het VAF bands: CN2 {1/2}, CN3 {1/3, 2/3}, CN4 {1/4..3/4}, CN>=5 broad;
    CN0/CN1 should show no hets (mismap floor). Per-bin contribution capped at
    ``max_sites_weight`` sites so SNP-dense bins don't overwhelm the depth channel.
    """
    from scipy.stats import beta as beta_dist

    emit = np.zeros((n_bins, max_cn + 1))
    het = (baf_vaf > het_lo) & (baf_vaf < het_hi)
    if not het.any():
        return emit
    pos, vaf = baf_pos[het], baf_vaf[het]
    bins = np.minimum(pos // bin_size, n_bins - 1)

    def band_loglik(v: np.ndarray, means: list[float]) -> np.ndarray:
        parts = [
            beta_dist.logpdf(
                np.clip(v, 1e-3, 1 - 1e-3), m * concentration, (1 - m) * concentration
            )
            for m in means
        ]
        stacked = np.stack(parts)
        m = stacked.max(axis=0)
        return m + np.log(np.exp(stacked - m).mean(axis=0))

    log_mismap = np.log(mismap)
    per_site = np.empty((len(vaf), max_cn + 1))
    for cn in range(max_cn + 1):
        if cn <= 1:
            per_site[:, cn] = log_mismap  # a clear het contradicts CN<=1
        elif cn == 2:
            per_site[:, cn] = band_loglik(vaf, [0.5])
        else:
            means = [i / cn for i in range(1, cn)]
            per_site[:, cn] = band_loglik(vaf, means)

    counts = np.zeros(n_bins)
    np.add.at(counts, bins, 1.0)
    weight = np.ones(len(vaf))
    dense = counts[bins] > max_sites_weight
    weight[dense] = max_sites_weight / counts[bins][dense]
    for cn in range(max_cn + 1):
        np.add.at(emit[:, cn], bins, per_site[:, cn] * weight)
    return emit


# ------------------------------------------------------------------------ segmentation


def segment_contig(
    ratio: np.ndarray,
    mask: np.ndarray,
    ploidy: int,
    max_cn: int,
    transition: float,
    baf_emit: np.ndarray | None = None,
    baf_weight: float = 1.0,
    bridge_max_gap_bins: int = 0,
    bridge_max_gap_frac: float = 0.05,
    bridge_max_ratio_diff: float = 0.15,
    extra_sd: np.ndarray | None = None,
) -> tuple[list[Segment], np.ndarray]:
    """Run the HMM on one contig; return non-baseline segments and the CN path.

    With ``bridge_max_gap_bins`` > 0, same-CN segments split by a short
    interruption are joined (see :func:`bridge_segments`); the path is the
    raw Viterbi path either way.
    """
    n = len(ratio)
    if n == 0:
        return [], np.array([], dtype=np.int32)
    sigma0 = fit_sigma0(ratio, mask, ploidy)
    emit = depth_log_emissions(ratio, mask, ploidy, max_cn, sigma0, extra_sd)
    if baf_emit is not None and baf_weight > 0:
        emit = emit + baf_weight * baf_emit

    k = max_cn + 1
    baseline = min(ploidy, max_cn)
    log_pi = np.full(k, np.log(1e-6))
    log_pi[baseline] = np.log(1 - 1e-6 * (k - 1))
    log_stay = float(np.log1p(-transition * (k - 1)))
    log_switch = float(np.log(transition))

    path = _viterbi(emit, log_stay, log_switch, log_pi)
    gamma = _forward_backward(emit, log_stay, log_switch, log_pi)

    segments: list[Segment] = []
    boundaries = np.flatnonzero(np.diff(path)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [n]))
    for s, e in zip(starts, ends, strict=True):
        cn = int(path[s])
        if cn == baseline:
            continue
        unmasked = ~mask[s:e]
        segments.append(
            Segment(
                start_bin=int(s),
                end_bin=int(e),
                cn=cn,
                posterior=_posterior(gamma[s:e, cn], unmasked),
                mean_ratio=float(ratio[s:e][unmasked].mean()) if unmasked.any() else float("nan"),
                n_bins=int(e - s),
                masked_frac=float(mask[s:e].mean()),
            )
        )
    n_raw = len(segments)
    segments = bridge_segments(
        segments, ratio, mask, gamma,
        max_gap_bins=bridge_max_gap_bins,
        max_gap_frac=bridge_max_gap_frac,
        max_ratio_diff=bridge_max_ratio_diff,
    )
    segments = [trim_masked_ends(seg, ratio, mask, gamma) for seg in segments]
    logger.debug(
        "HMM: %d bins, sigma0=%.3f, %d non-baseline segments (%d after bridging)",
        n, sigma0, n_raw, len(segments),
    )
    return segments, path


def _posterior(g: np.ndarray, informative: np.ndarray) -> float:
    """Mean posterior of the called state over the informative (unmasked) bins.

    A masked bin has a flat emission: its posterior only echoes its distance
    to the nearest informative bins, and a 150 kb segmental-duplication run
    inside a 16p11.2 deletion dragged the mean under ``min_cn_posterior``."""
    if informative.any():
        return float(g[informative].mean())
    return float(g.mean())


def trim_masked_ends(seg: Segment, ratio: np.ndarray, mask: np.ndarray,
                     gamma: np.ndarray) -> Segment:
    """Drop the masked bins at either end of a segment.

    The HMM carries a CN state through masked (uninformative) bins, so a
    deletion flanked by segmental duplications (16p11.2, 15q11.2, 1q21.1)
    ran on through them and its call inherited the duplications' MASKFRAC and
    SEGDUP_OVL (LowMap, MaxN). The call keeps its informative core; the
    trimmed bins widen CIPOS/CIEND. A segment with no unmasked bin is
    returned unchanged (MaxN)."""
    s, e = seg.start_bin, seg.end_bin
    while s < e and mask[s]:
        s += 1
    while e > s and mask[e - 1]:
        e -= 1
    if s == e or (s == seg.start_bin and e == seg.end_bin):
        return seg
    unmasked = ~mask[s:e]
    return Segment(
        start_bin=s, end_bin=e, cn=seg.cn,
        posterior=_posterior(gamma[s:e, seg.cn], unmasked),
        mean_ratio=float(ratio[s:e][unmasked].mean()),
        n_bins=e - s,
        masked_frac=float(mask[s:e].mean()),
        trim_left=seg.trim_left + s - seg.start_bin,
        trim_right=seg.trim_right + seg.end_bin - e,
    )


def bridge_segments(
    segments: list[Segment],
    ratio: np.ndarray,
    mask: np.ndarray,
    gamma: np.ndarray,
    max_gap_bins: int,
    max_gap_frac: float,
    max_ratio_diff: float,
) -> list[Segment]:
    """Join same-CN segments separated by a short interruption.

    The Viterbi path leaves a CN1 run for a single outlier bin (ratio ~1.3,
    cheaper as CN3 than as a 12-sigma CN1 bin) or for the 2 kb CN0 of a
    common deletion on the remaining allele, so one 6 Mb terminal deletion
    came out as 8 abutting calls. Segments A and B (same CN, in path order)
    are joined when the bins between them number at most ``max_gap_bins``
    and at most ``max_gap_frac`` of the shorter of the two sides, and their
    mean ratios agree within ``max_ratio_diff``. Joins chain (A+B, then
    +C). The joined segment spans the interruption; its posterior and mean
    ratio are taken over the flanking bins only. Segments inside the
    interruption are returned unchanged, nested in the joined one.
    """
    if max_gap_bins <= 0 or len(segments) < 2:
        return segments
    out: list[Segment] = []
    absorbed = [False] * len(segments)
    for i, seg in enumerate(segments):
        if absorbed[i]:
            continue
        run = [seg]
        run_bins = seg.n_bins
        run_ratio = seg.mean_ratio
        j = i + 1
        while j < len(segments):
            nxt = segments[j]
            gap = nxt.start_bin - run[-1].end_bin
            if gap > max_gap_bins:
                break
            if nxt.cn != seg.cn:
                j += 1  # a segment inside the interruption (or past it: gap check above)
                continue
            ok = (
                gap > 0
                and gap <= max_gap_frac * min(run_bins, nxt.n_bins)
                and np.isfinite(run_ratio)
                and np.isfinite(nxt.mean_ratio)
                and abs(run_ratio - nxt.mean_ratio) <= max_ratio_diff
            )
            if not ok:
                break
            run_ratio = (run_ratio * run_bins + nxt.mean_ratio * nxt.n_bins) / (
                run_bins + nxt.n_bins
            )
            run_bins += nxt.n_bins
            run.append(nxt)
            absorbed[j] = True
            j += 1
        if len(run) == 1:
            out.append(seg)
            continue
        start, end = run[0].start_bin, run[-1].end_bin
        flank = np.concatenate([np.arange(s.start_bin, s.end_bin) for s in run])
        unmasked = flank[~mask[flank]]
        out.append(
            Segment(
                start_bin=start,
                end_bin=end,
                cn=seg.cn,
                posterior=_posterior(gamma[flank, seg.cn], ~mask[flank]),
                mean_ratio=float(ratio[unmasked].mean()) if len(unmasked) else float("nan"),
                n_bins=end - start,
                masked_frac=float(mask[start:end].mean()),
            )
        )
    out.sort(key=lambda s: (s.start_bin, s.end_bin))
    return out
