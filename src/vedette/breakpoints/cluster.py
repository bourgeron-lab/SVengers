"""Clustering of discordant-pair and split-read evidence into junctions.

A *junction* is a putative adjacency (two breakpoints on the same contig for
DEL/tandem-DUP). PE clusters give imprecise intervals; SR pileups refine them
to base pairs. SR-only junctions are also formed (small events below insert
resolution) but flagged as such for the downstream filter.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from ..config import BreakpointConfig
from ..extract.store import PE_DEL, PE_DUP, SR_LEFT_CLIP, SR_RIGHT_CLIP, InsertModel

logger = logging.getLogger(__name__)


@dataclass
class PECluster:
    orient: int
    n: int
    mean_mapq: float
    pos1_lo: int
    pos1_hi: int
    pos2_lo: int
    pos2_hi: int


@dataclass
class SRBreakpoint:
    pos: int
    side: int          # SR_LEFT_CLIP / SR_RIGHT_CLIP
    n: int
    mean_mapq: float
    max_cliplen: int


@dataclass
class Junction:
    svtype: str        # "DEL" | "DUP"
    start: int         # 0-based breakpoint (left)
    end: int           # 0-based breakpoint (right, exclusive)
    cipos: tuple[int, int]
    ciend: tuple[int, int]
    pe: int
    sr_start: int
    sr_end: int
    precise: bool
    homlen: int = -1   # reference micro-homology at the junction (-1 = not measured)

    @property
    def sr(self) -> int:
        return self.sr_start + self.sr_end

    @property
    def size(self) -> int:
        return self.end - self.start


def cluster_pe(
    pe: dict[str, np.ndarray],
    insert: InsertModel,
    config: BreakpointConfig,
    orient_wanted: int,
) -> list[PECluster]:
    """Greedy single-linkage clustering of same-orientation pairs sorted by pos1.

    A pair joins the open cluster if its pos1 gap is within the insert spread
    AND its pos2 is compatible with the cluster's running pos2 range.
    """
    sel = pe["orient"] == orient_wanted
    if not sel.any():
        return []
    pos1 = pe["pos1"][sel].astype(np.int64)
    pos2 = pe["pos2"][sel].astype(np.int64)
    mapq = pe["mapq"][sel].astype(np.float64)
    # already sorted by pos1 in the store; enforce for safety
    order = np.argsort(pos1, kind="stable")
    pos1, pos2, mapq = pos1[order], pos2[order], mapq[order]

    gap = insert.mu + 3 * insert.sigma
    clusters: list[PECluster] = []
    start_idx = 0
    p2_lo = p2_hi = pos2[0]

    def _flush(a: int, b: int, lo: int, hi: int) -> None:
        n = b - a
        if n >= config.min_pe:
            clusters.append(
                PECluster(
                    orient=orient_wanted,
                    n=n,
                    mean_mapq=float(mapq[a:b].mean()),
                    pos1_lo=int(pos1[a]),
                    pos1_hi=int(pos1[b - 1]),
                    pos2_lo=int(lo),
                    pos2_hi=int(hi),
                )
            )

    for i in range(1, len(pos1)):
        compatible = (
            pos1[i] - pos1[i - 1] <= gap
            and min(abs(pos2[i] - p2_lo), abs(pos2[i] - p2_hi)) <= gap
        )
        if compatible:
            p2_lo, p2_hi = min(p2_lo, pos2[i]), max(p2_hi, pos2[i])
        else:
            _flush(start_idx, i, p2_lo, p2_hi)
            start_idx = i
            p2_lo = p2_hi = pos2[i]
    _flush(start_idx, len(pos1), p2_lo, p2_hi)
    return clusters


def cluster_sr(
    sr: dict[str, np.ndarray], config: BreakpointConfig
) -> list[SRBreakpoint]:
    """Pile up clip positions (per side) within the position tolerance."""
    out: list[SRBreakpoint] = []
    for side in (SR_LEFT_CLIP, SR_RIGHT_CLIP):
        sel = sr["side"] == side
        if not sel.any():
            continue
        pos = sr["pos"][sel].astype(np.int64)
        mapq = sr["mapq"][sel].astype(np.float64)
        clip = sr["cliplen"][sel].astype(np.int64)
        order = np.argsort(pos, kind="stable")
        pos, mapq, clip = pos[order], mapq[order], clip[order]

        start = 0
        for i in range(1, len(pos) + 1):
            if i == len(pos) or pos[i] - pos[i - 1] > config.sr_pos_tolerance:
                n = i - start
                if n >= config.min_sr:
                    out.append(
                        SRBreakpoint(
                            pos=int(np.median(pos[start:i])),
                            side=side,
                            n=n,
                            mean_mapq=float(mapq[start:i].mean()),
                            max_cliplen=int(clip[start:i].max()),
                        )
                    )
                start = i
    out.sort(key=lambda b: b.pos)
    return out


def _nearest_sr(
    breakpoints: list[SRBreakpoint], side: int, lo: int, hi: int
) -> SRBreakpoint | None:
    """Best SR pileup of the wanted side within [lo, hi] (most reads wins)."""
    best = None
    for bp in breakpoints:
        if bp.side != side or not (lo <= bp.pos <= hi):
            continue
        if best is None or bp.n > best.n:
            best = bp
    return best


def make_junctions(
    pe: dict[str, np.ndarray],
    sr: dict[str, np.ndarray],
    insert: InsertModel,
    config: BreakpointConfig,
) -> list[Junction]:
    """Combine PE clusters and SR pileups into DEL/DUP junction candidates.

    Geometry (leftmost read recorded, FR library):
    - DEL: +/- pair with stretched insert. Left breakpoint lies right of the
      forward reads (pos1_hi .. pos1_hi + mu); right breakpoint left of the
      reverse reads (pos2_lo - mu .. pos2_lo + read length). At the left
      breakpoint reads are clipped on their RIGHT; at the right breakpoint on
      their LEFT.
    - Tandem DUP: everted -/+ pair. The junction joins the DUP's END (right
      clip side, near pos2) back to its START (left clip side, near pos1):
      start breakpoint in [pos1_lo - mu, pos1_hi], end in [pos2_lo, pos2_hi + mu].
    """
    sr_bps = cluster_sr(sr, config)
    junctions: list[Junction] = []
    mu = int(insert.mu) or 500

    for orient, svtype in ((PE_DEL, "DEL"), (PE_DUP, "DUP")):
        for cl in cluster_pe(pe, insert, config, orient):
            if orient == PE_DEL:
                start_lo, start_hi = cl.pos1_hi, cl.pos1_hi + mu
                end_lo, end_hi = max(cl.pos2_lo - mu, 0), cl.pos2_lo + 150
            else:
                start_lo, start_hi = max(cl.pos1_lo - mu, 0), cl.pos1_hi + 150
                end_lo, end_hi = cl.pos2_lo, cl.pos2_hi + mu

            sr_at_start = _nearest_sr(sr_bps, SR_RIGHT_CLIP, start_lo, start_hi)
            sr_at_end = _nearest_sr(sr_bps, SR_LEFT_CLIP, end_lo, end_hi)

            start = sr_at_start.pos if sr_at_start else (start_lo + start_hi) // 2
            end = sr_at_end.pos if sr_at_end else (end_lo + end_hi) // 2
            if end - start < config.min_size or end - start > config.max_junction_size:
                continue
            precise = sr_at_start is not None and sr_at_end is not None
            junctions.append(
                Junction(
                    svtype=svtype,
                    start=start,
                    end=end,
                    cipos=(0, 0) if sr_at_start else (start_lo - start, start_hi - start),
                    ciend=(0, 0) if sr_at_end else (end_lo - end, end_hi - end),
                    pe=cl.n,
                    sr_start=sr_at_start.n if sr_at_start else 0,
                    sr_end=sr_at_end.n if sr_at_end else 0,
                    precise=precise,
                )
            )

    junctions.sort(key=lambda j: (j.start, j.end))
    logger.debug("junctions: %d (from %d SR pileups)", len(junctions), len(sr_bps))
    return junctions
