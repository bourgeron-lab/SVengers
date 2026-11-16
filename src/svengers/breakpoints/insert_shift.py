"""Sub-insert deletion detection from the local insert-size shift.

Proper pairs spanning a deletion smaller than the insert size are individually
concordant (their stretch is below the discordance threshold) but collectively
shifted: the local mean template length rises by the deleted length. This is
the SurVIndel2 insight, and it is the only usable signal for the 50-300 bp
deletions that neither clip cleanly (VNTR contractions scatter the alignment
into micro-indels) nor move a 1 kb depth bin.

Detection: per-bin z-score of (local mean insert - mu) against sigma/sqrt(n),
runs of consecutive shifted bins -> one candidate. The weighted mean shift d
estimates the HAPLOID deleted length times the carrier fraction; assuming the
common heterozygous case, size ~ 2d. Positions are imprecise by construction
(the shifted region spans ~ one insert length): candidates carry wide CIs and
SOURCE=INSERT, and are pre-genotyped 0/1 with a z-derived GQ — the calibrated
filter decides their fate from ins_z / ins_pairs features.
"""

from __future__ import annotations

import logging

import numpy as np

from ..config import BreakpointConfig
from ..extract.store import InsertModel
from ..variants import SOURCE_INSERT, CNVCall

logger = logging.getLogger(__name__)


def detect_insert_shift_dels(
    contig: str,
    ins_sum: np.ndarray,
    ins_n: np.ndarray,
    bin_size: int,
    insert: InsertModel,
    config: BreakpointConfig,
    min_n_per_bin: int = 3,
) -> list[CNVCall]:
    if insert.sigma <= 0 or not len(ins_sum):
        return []
    n = ins_n.astype(np.float64)
    ok = n >= min_n_per_bin
    mean = np.where(ok, ins_sum / np.maximum(n, 1), insert.mu)
    delta = mean - insert.mu
    z = np.where(ok, delta / (insert.sigma / np.sqrt(np.maximum(n, 1))), 0.0)

    min_shift = config.min_size / 2.0 * 0.8  # het 50bp DEL shifts by ~25
    shifted = (z >= config.insert_shift_z) & (delta >= min_shift)

    calls: list[CNVCall] = []
    idx = np.flatnonzero(shifted)
    if not len(idx):
        return calls
    # group consecutive shifted bins (allowing single-bin dropouts)
    run_start = idx[0]
    prev = idx[0]
    for b in list(idx[1:]) + [-1]:
        if b != -1 and b - prev <= 2:
            prev = b
            continue
        run_bins = np.arange(run_start, prev + 1)
        run_bins = run_bins[shifted[run_bins]]
        if len(run_bins) >= config.insert_shift_min_bins:
            n_pairs = int(n[run_bins].sum())
            if n_pairs >= config.insert_shift_min_pairs:
                d = float((delta[run_bins] * n[run_bins]).sum() / n[run_bins].sum())
                size_est = max(int(round(2 * d)), config.min_size)
                center = int((run_bins[0] + run_bins[-1] + 1) / 2 * bin_size)
                start = max(center - size_est // 2, 0)
                zmax = float(z[run_bins].max())
                ci = int(insert.mu)
                call = CNVCall(
                    contig=contig,
                    start=start,
                    end=start + size_est,
                    svtype="DEL",
                    source=SOURCE_INSERT,
                    cipos=(-ci, ci),
                    ciend=(-ci, ci),
                    precise=False,
                    ins_z=round(zmax, 2),
                    ins_pairs=n_pairs,
                    gt=(0, 1),
                    gq=int(min(60, 3 * zmax)),
                )
                calls.append(call)
        if b == -1:
            break
        run_start = prev = b
    logger.debug("%s: %d insert-shift DEL candidates", contig, len(calls))
    return calls


def dedupe_against(
    calls: list[CNVCall], existing: list[CNVCall], min_overlap: float = 0.3
) -> list[CNVCall]:
    """Drop insert-shift candidates already covered by junction/depth calls."""
    kept = []
    for c in calls:
        redundant = False
        for e in existing:
            if e.svtype != "DEL" or e.contig != c.contig:
                continue
            inter = min(c.end, e.end) - max(c.start, e.start)
            if inter > 0 and inter >= min_overlap * min(c.size, e.size):
                redundant = True
                break
        if not redundant:
            kept.append(c)
    return kept
