"""Sample-level QC metrics computed at calling time.

Wave score: some libraries show megabase-scale sinusoidal coverage
oscillations ("waves", insert-size/GC batch artifacts) that survive GC
correction. Wavy samples fragment the CN-HMM into spurious segments, flood the
hidden-split realignment with candidate windows (5x median calling runtime on
GHFC) and carry a higher false-call rate — the score is stored in the evidence
store, exposed in the cohort report and fed to the filter model as a covariate.

The metric follows the GHFC wave-detector's dominant component: the fraction
of spectral energy at low frequencies. Depth ratio is cut into 10 Mb windows;
per window, the linearly-detrended signal's FFT energy at periods >= 200 kb is
divided by the total energy; per-contig score is the median over windows and
the sample score the mean over contigs. Range [0, 1]; flat noise ~= the
fraction of spectrum bins in the wave band (a few percent), wavy GHFC samples
reach 0.5+ on the corrected ratio.
"""

from __future__ import annotations

import numpy as np

WINDOW_BP = 10_000_000   # spectral window
# Report/QC warning threshold. Calibrated on 592 GHFC genomes against the GHFC
# wave-detector labels (AUC 0.946): 0.145 = p95 of normal samples — mirrors the
# detector's own 95th-percentile calibration — sensitivity 0.80 on labeled
# waves. On the CORRECTED ratio the dynamic range is compressed (GC correction
# absorbs part of the artifact), hence thresholds well below the raw-coverage
# scores. See docs/wave_qc_validation.md.
WAVE_WARN = 0.145
WAVE_MAX_K = 50          # FFT bins 1..50 = periods >= WINDOW_BP/50 = 200 kb
MIN_VALID_FRAC = 0.5     # a window needs this fraction of unmasked bins


def contig_wave_score(
    ratio: np.ndarray, mask: np.ndarray, bin_size: int
) -> float | None:
    """Low-frequency spectral energy fraction of one contig's copy-ratio track.

    Masked bins are filled with the window mean (spectrally neutral) rather
    than interpolated, so long masked stretches cannot fabricate waves.
    """
    n_win = WINDOW_BP // bin_size
    if n_win < 4 * WAVE_MAX_K or len(ratio) < n_win:
        return None
    scores = []
    x = np.arange(n_win, dtype=np.float64)
    for lo in range(0, len(ratio) - n_win + 1, n_win // 2):
        win = ratio[lo:lo + n_win].astype(np.float64)
        ok = ~mask[lo:lo + n_win] & np.isfinite(win)
        if ok.mean() < MIN_VALID_FRAC:
            continue
        w = win.copy()
        w[~ok] = win[ok].mean()
        if w.std() < 1e-6:
            continue
        # linear detrend (numpy-only)
        slope, intercept = np.polyfit(x, w, 1)
        w -= slope * x + intercept
        power = np.abs(np.fft.rfft(w)) ** 2
        total = float(power[1:].sum())
        if total <= 0:
            continue
        scores.append(float(power[1:WAVE_MAX_K + 1].sum()) / total)
    return float(np.median(scores)) if scores else None


def sample_wave_score(
    ratio_by_contig: dict[str, np.ndarray],
    mask_by_contig: dict[str, np.ndarray],
    bin_size: int,
    min_contigs: int = 10,
) -> tuple[float | None, dict[str, float]]:
    """(genome score, per-contig scores). None when too few contigs qualify
    (e.g. single-contig smoke runs — the per-contig dict is still returned)."""
    per_contig = {}
    for name, ratio in ratio_by_contig.items():
        s = contig_wave_score(ratio, mask_by_contig[name], bin_size)
        if s is not None:
            per_contig[name] = s
    if len(per_contig) < min_contigs:
        return None, per_contig
    return float(np.mean(list(per_contig.values()))), per_contig
