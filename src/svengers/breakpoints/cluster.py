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
    cig: int = 0       # reads carrying the event as a CIGAR deletion

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
    """Pile up clip positions (per side) within the position tolerance.

    Vectorised over the groups (a whole contig's clips are millions; the
    per-clip loop was 70 % of `svengers bnd`): same pileups, positions and
    order as the loop kept in the tests as reference."""
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

        starts = np.concatenate([[0], np.flatnonzero(np.diff(pos) > config.sr_pos_tolerance) + 1])
        n = np.diff(np.concatenate([starts, [len(pos)]]))
        keep = n >= config.min_sr
        mapq_sum = np.add.reduceat(mapq, starts)[keep]
        clip_max = np.maximum.reduceat(clip, starts)[keep]
        a, n = starts[keep], n[keep]
        # the median of a sorted group: its middle clip, or the mean of the two
        med = (pos[a + (n - 1) // 2] + pos[a + n // 2]) / 2
        out += [
            SRBreakpoint(pos=int(m), side=side, n=int(k), mean_mapq=float(q / k),
                         max_cliplen=int(c))
            for m, k, q, c in zip(med, n, mapq_sum, clip_max, strict=True)
        ]
    out.sort(key=lambda b: b.pos)
    return out


def cluster_cigar_dels(
    cig: dict[str, np.ndarray], config: BreakpointConfig, pos_tol: int = 10
) -> list[Junction]:
    """Cluster large CIGAR deletions into precise DEL junctions.

    Aligners left-align indels, so reads carrying the same deletion agree on
    (pos, len) almost exactly; a small tolerance absorbs tandem-repeat jitter.
    """
    if not len(cig.get("pos", ())):
        return []
    pos = cig["pos"].astype(np.int64)
    svlen = cig["svlen"].astype(np.int64)
    order = np.lexsort((svlen, pos))
    pos, svlen = pos[order], svlen[order]

    junctions: list[Junction] = []
    start_idx = 0
    for i in range(1, len(pos) + 1):
        boundary = (
            i == len(pos)
            or pos[i] - pos[start_idx] > pos_tol
            or abs(svlen[i] - svlen[start_idx]) > max(pos_tol, 0.1 * svlen[start_idx])
        )
        if not boundary:
            continue
        n = i - start_idx
        if n >= config.min_cigar_support:
            p = int(np.median(pos[start_idx:i]))
            ln = int(np.median(svlen[start_idx:i]))
            if config.min_size <= ln <= config.max_junction_size:
                junctions.append(
                    Junction(
                        svtype="DEL", start=p, end=p + ln,
                        cipos=(0, 0), ciend=(0, 0),
                        pe=0, sr_start=0, sr_end=0,
                        precise=True, cig=n,
                    )
                )
        start_idx = i
    return junctions


def _nearest_sr(
    breakpoints: list[SRBreakpoint], side: int, lo: int, hi: int
) -> SRBreakpoint | None:
    """Best SR pileup of the wanted side within [lo, hi] (most reads wins;
    the first of equals in position order). Reference for _SRIndex."""
    best = None
    for bp in breakpoints:
        if bp.side != side or not (lo <= bp.pos <= hi):
            continue
        if best is None or bp.n > best.n:
            best = bp
    return best


class _SRIndex:
    """_nearest_sr by binary search on per-side, position-sorted pileups.

    The linear scan ran once per PE cluster over every pileup of the contig:
    quadratic, and a noisy genome has millions of pileups (a batch of five
    noisy dataset_paper genomes: 1.9 M on chr20, 7 of the 13 min of its chr20
    call, 5-9 h per genome). Same answer: the pileups of one side keep their relative
    order, so the first maximum in the window is the same pileup.
    """

    def __init__(self, breakpoints: list[SRBreakpoint]):
        self._by_side: dict[int, tuple[np.ndarray, list[SRBreakpoint]]] = {}
        for side in {bp.side for bp in breakpoints}:
            bps = [bp for bp in breakpoints if bp.side == side]
            order = sorted(range(len(bps)), key=lambda i: bps[i].pos)  # stable
            bps = [bps[i] for i in order]
            self._by_side[side] = (np.array([bp.pos for bp in bps], dtype=np.int64), bps)

    def nearest(self, side: int, lo: int, hi: int) -> SRBreakpoint | None:
        entry = self._by_side.get(side)
        if entry is None:
            return None
        pos, bps = entry
        i = int(np.searchsorted(pos, lo, "left"))
        j = int(np.searchsorted(pos, hi, "right"))
        best = None
        for bp in bps[i:j]:
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
    sr_index = _SRIndex(sr_bps)
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

            sr_at_start = sr_index.nearest(SR_RIGHT_CLIP, start_lo, start_hi)
            sr_at_end = sr_index.nearest(SR_LEFT_CLIP, end_lo, end_hi)

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
