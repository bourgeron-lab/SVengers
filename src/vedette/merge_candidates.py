"""Reconciliation of HMM depth segments and PE/SR junctions into candidates.

- depth segment + concordant junction  -> one call, junction coordinates
  (base-pair precise when both ends have split reads), source=BOTH;
- junction alone -> call kept below the size where depth is underpowered,
  source=JUNCTION;
- depth segment alone -> imprecise call, source=DEPTH (large CNVs without
  clean breakpoints, e.g. in repeats).
"""

from __future__ import annotations

import logging

import numpy as np

from .breakpoints.cluster import Junction
from .config import Config
from .segment.hmm import Segment
from .variants import SOURCE_BOTH, SOURCE_DEPTH, SOURCE_JUNCTION, CNVCall

logger = logging.getLogger(__name__)


def segment_to_call(
    contig: str, seg: Segment, bin_size: int, ploidy: int, config: Config
) -> CNVCall:
    svtype = "DEL" if seg.cn < ploidy else "DUP"
    call = CNVCall(
        contig=contig,
        start=seg.start_bin * bin_size,
        end=seg.end_bin * bin_size,
        svtype=svtype,
        source=SOURCE_DEPTH,
        cn=seg.cn,
        cn_posterior=seg.posterior,
        mean_ratio=seg.mean_ratio,
        n_bins=seg.n_bins,
        masked_frac=seg.masked_frac,
        cipos=(-bin_size, bin_size),
        ciend=(-bin_size, bin_size),
        precise=False,
    )
    if seg.n_bins < config.segment.min_bins:
        call.filters.append("NoEvidence")
    if seg.masked_frac > config.segment.max_masked_frac:
        call.filters.append("MaxN")
    if seg.posterior < config.segment.min_cn_posterior:
        call.filters.append("LowPosterior")
    if call.size < config.breakpoint.min_size:
        call.filters.append("MinSize")
    return call


def _boundary_tolerance(size: int, bin_size: int) -> int:
    return max(2 * bin_size, int(0.1 * size))


def reconcile(
    contig: str,
    segments: list[Segment],
    junctions: list[Junction],
    bin_size: int,
    ploidy: int,
    config: Config,
) -> list[CNVCall]:
    calls: list[CNVCall] = []
    seg_calls = [
        segment_to_call(contig, seg, bin_size, ploidy, config)
        for seg in segments
        if seg.cn != ploidy
    ]
    used_junction = np.zeros(len(junctions), dtype=bool)

    for call in seg_calls:
        tol = _boundary_tolerance(call.size, bin_size)
        best_idx, best_dist = -1, tol + 1
        for j, junc in enumerate(junctions):
            if used_junction[j] or junc.svtype != call.svtype:
                continue
            dist = max(abs(junc.start - call.start), abs(junc.end - call.end))
            if dist < best_dist:
                best_idx, best_dist = j, dist
        if best_idx >= 0:
            junc = junctions[best_idx]
            used_junction[best_idx] = True
            call.start, call.end = junc.start, junc.end
            call.cipos, call.ciend = junc.cipos, junc.ciend
            call.pe_support = junc.pe
            call.sr_support = junc.sr
            call.cig_support = junc.cig
            call.precise = junc.precise
            call.source = SOURCE_BOTH
            if junc.homlen >= 0:
                call.info_extra["HOMLEN"] = junc.homlen
            # junction-confirmed calls are exempt from depth-only caveats
            call.filters = [
                f for f in call.filters if f not in ("LowPosterior", "NoEvidence")
            ]
        calls.append(call)

    max_no_depth = config.breakpoint.junction_only_max_size
    for j, junc in enumerate(junctions):
        if used_junction[j] or junc.size > max_no_depth:
            continue
        call = CNVCall(
            contig=contig,
            start=junc.start,
            end=junc.end,
            svtype=junc.svtype,
            source=SOURCE_JUNCTION,
            cipos=junc.cipos,
            ciend=junc.ciend,
            pe_support=junc.pe,
            sr_support=junc.sr,
            cig_support=junc.cig,
            precise=junc.precise,
        )
        if junc.homlen >= 0:
            call.info_extra["HOMLEN"] = junc.homlen
        if junc.size < config.breakpoint.min_size:
            call.filters.append("MinSize")
        calls.append(call)

    calls.sort(key=lambda c: (c.start, c.end))
    n_both = sum(1 for c in calls if c.source == SOURCE_BOTH)
    logger.debug(
        "%s: %d calls (%d depth+junction, %d junction-only)",
        contig, len(calls), n_both,
        sum(1 for c in calls if c.source == SOURCE_JUNCTION),
    )
    return calls
