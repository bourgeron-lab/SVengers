"""Hidden-split-read recovery: local realignment of soft-clip sequences.

The SurVIndel2 insight, reimplemented on our evidence: a soft-clipped read is
half a split read whose other half the aligner dropped. Realigning the clipped
sequence against the local reference (± a window) turns the clip into a
base-pair-precise junction — and the DIRECTION of the hit gives the SV type:

    RIGHT clip at p, clip start realigns at q:   q > p  -> DEL(p, q)
                                                 q < p  -> tandem DUP(q, p)
    LEFT  clip at p, clip end   realigns at q:   q < p  -> DEL(q, p)
                                                 q > p  -> tandem DUP(p, q)

This recovers tandem duplications SMALLER THAN THE INSERT SIZE, which produce
no everted read pairs and are the historic blind spot of short-read callers
(all tools <=0.47 DUP sensitivity, worst at 1-10 kb).

Alignment is a fit alignment (query-global, reference-local) in numba with
O(reference) memory and start-position tracking, so no external aligner
dependency. Junction micro-homology is measured on the reference and reported
(HOMLEN): high homology marks tandem-repeat junctions for the filter.

Reverse-complement hits (inversions) are deliberately out of scope until v2.
"""

from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np

from ..config import BreakpointConfig
from ..extract.store import SR_RIGHT_CLIP
from .cluster import Junction, cluster_sr

logger = logging.getLogger(__name__)

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        def wrap(f):
            return f
        return wrap if not (args and callable(args[0])) else args[0]

_MATCH, _MISMATCH, _GAP = 2, -3, -4


@njit(cache=True)
def _fit_align(query: np.ndarray, ref: np.ndarray):  # pragma: no cover - numba
    """Fit alignment: full query against any substring of ref.

    Returns (score, ref_start, ref_end) of the best placement (0-based,
    half-open on ref). O(len(ref)) memory; alignment start tracked per cell.
    """
    m, n = len(query), len(ref)
    prev = np.zeros(n + 1, dtype=np.int32)
    prev_start = np.arange(n + 1, dtype=np.int32)  # alignment begins at ref j
    cur = np.zeros(n + 1, dtype=np.int32)
    cur_start = np.zeros(n + 1, dtype=np.int32)

    for i in range(1, m + 1):
        cur[0] = prev[0] + _GAP
        cur_start[0] = 0
        qc = query[i - 1]
        for j in range(1, n + 1):
            diag = prev[j - 1] + (_MATCH if ref[j - 1] == qc else _MISMATCH)
            up = prev[j] + _GAP      # consume query, gap in ref
            left = cur[j - 1] + _GAP  # consume ref, gap in query
            if diag >= up and diag >= left:
                cur[j] = diag
                cur_start[j] = prev_start[j - 1]
            elif up >= left:
                cur[j] = up
                cur_start[j] = prev_start[j]
            else:
                cur[j] = left
                cur_start[j] = cur_start[j - 1]
        prev, cur = cur, prev
        prev_start, cur_start = cur_start, prev_start

    best_j, best_score = 0, np.int32(-2147483647)
    for j in range(n + 1):
        if prev[j] > best_score:
            best_score = prev[j]
            best_j = j
    return int(best_score), int(prev_start[best_j]), int(best_j)


def _seq_to_arr(seq: bytes) -> np.ndarray:
    return np.frombuffer(seq.upper(), dtype=np.uint8)


def microhomology(ref: np.ndarray, start: int, end: int, cap: int = 50) -> int:
    """Reference homology across the adjacency (start, end): forward + backward."""
    n = len(ref)
    fwd = 0
    while (
        fwd < cap and start + fwd < n and end + fwd < n
        and ref[start + fwd] == ref[end + fwd]
    ):
        fwd += 1
    bwd = 0
    while (
        bwd < cap and start - 1 - bwd >= 0 and end - 1 - bwd >= 0
        and ref[start - 1 - bwd] == ref[end - 1 - bwd]
    ):
        bwd += 1
    return fwd + bwd


def recover_junctions(
    sr: dict[str, np.ndarray],
    sr_seq: list[bytes],
    ref: np.ndarray,
    config: BreakpointConfig,
) -> list[Junction]:
    """Realign clip sequences at SR pileups -> precise DEL/DUP junctions."""
    if not len(sr.get("pos", ())) or not sr_seq:
        return []
    pileups = cluster_sr(sr, config)
    if not pileups:
        return []
    # runtime guard: strongest pileups first, stop at the alignment budget
    pileups = sorted(pileups, key=lambda b: (-b.n, -b.max_cliplen))
    budget = config.resolve_max_alignments
    n_aligned = 0
    n_skipped_pileups = 0

    pos_all = sr["pos"].astype(np.int64)
    side_all = sr["side"]
    clip_all = sr["cliplen"].astype(np.int64)
    window = config.resolve_window

    # position-sorted clip index per side: pileup membership by binary search,
    # not a full-array scan per pileup (quadratic on junction-dense genomes —
    # profiled at hours/contig on wavy GHFC samples)
    by_side: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for side in np.unique(side_all):
        idx = np.flatnonzero(
            (side_all == side) & (clip_all >= config.resolve_min_clip)
        )
        idx = idx[np.argsort(pos_all[idx], kind="stable")]
        by_side[int(side)] = (pos_all[idx], idx)

    # candidate adjacencies: (svtype, start, end) -> [homlen-capable evidence count]
    votes: dict[tuple[str, int, int], int] = defaultdict(int)

    for bp in pileups:
        if n_aligned >= budget:
            n_skipped_pileups += 1
            continue
        # clips belonging to this pileup, longest first, capped
        if bp.side not in by_side:
            continue
        sorted_pos, sorted_idx = by_side[bp.side]
        lo = np.searchsorted(sorted_pos, bp.pos - config.sr_pos_tolerance, "left")
        hi = np.searchsorted(sorted_pos, bp.pos + config.sr_pos_tolerance, "right")
        sel = sorted_idx[lo:hi]
        if len(sel) == 0:
            continue
        sel = sel[np.argsort(-clip_all[sel])][: config.resolve_max_clips]
        w_lo = max(bp.pos - window, 0)
        w_hi = min(bp.pos + window, len(ref))
        ref_win = ref[w_lo:w_hi]
        if len(ref_win) < 50:
            continue

        for idx in sel:
            seq = sr_seq[idx]
            if len(seq) < config.resolve_min_clip:
                continue
            query = _seq_to_arr(seq)
            n_aligned += 1
            score, rs, re_ = _fit_align(query, ref_win)
            if score < _MATCH * config.resolve_min_identity * len(query):
                continue
            p = int(pos_all[idx])
            if bp.side == SR_RIGHT_CLIP:
                q = w_lo + rs      # where the post-junction sequence starts
                if q > p + config.min_size:
                    key = ("DEL", p, q)
                elif q < p - config.min_size:
                    key = ("DUP", q, p)
                else:
                    continue
            else:
                q = w_lo + re_     # where the pre-junction sequence ends
                if q < p - config.min_size:
                    key = ("DEL", q, p)
                elif q > p + config.min_size:
                    key = ("DUP", p, q)
                else:
                    continue
            votes[key] += 1

    # cluster votes within a small tolerance and keep supported junctions
    junctions: list[Junction] = []
    by_type: dict[str, list[tuple[int, int, int]]] = {"DEL": [], "DUP": []}
    for (svtype, start, end), n in votes.items():
        by_type[svtype].append((start, end, n))
    tol = 15
    for svtype, cands in by_type.items():
        cands.sort()
        used = [False] * len(cands)
        for i, (s, e, _n) in enumerate(cands):
            if used[i]:
                continue
            support, s_list, e_list = 0, [], []
            for k in range(i, len(cands)):
                s2, e2, n2 = cands[k]
                if s2 - s > tol:
                    break
                if not used[k] and abs(e2 - e) <= tol:
                    used[k] = True
                    support += n2
                    s_list += [s2] * n2
                    e_list += [e2] * n2
            if support < config.resolve_min_support:
                continue
            start_c = int(np.median(s_list))
            end_c = int(np.median(e_list))
            junctions.append(
                Junction(
                    svtype=svtype,
                    start=start_c,
                    end=end_c,
                    cipos=(0, 0),
                    ciend=(0, 0),
                    pe=0,
                    sr_start=support,
                    sr_end=support,
                    precise=True,
                    homlen=microhomology(ref, start_c, end_c),
                )
            )
    junctions.sort(key=lambda j: (j.start, j.end))
    if n_skipped_pileups:
        logger.warning(
            "hidden splits: alignment budget %d exhausted, %d/%d weakest "
            "pileups skipped (wavy/noisy depth profile?)",
            budget, n_skipped_pileups, len(pileups),
        )
    logger.debug(
        "hidden splits: %d pileups -> %d junctions (%d DUP)",
        len(pileups), len(junctions),
        sum(1 for j in junctions if j.svtype == "DUP"),
    )
    return junctions


def merge_junction_sets(
    primary: list[Junction], recovered: list[Junction], tol: int = 50
) -> list[Junction]:
    """Fold realignment-recovered junctions into the PE/SR junction set.

    A recovered junction matching an existing one refines its coordinates
    (realignment is base-precise) and adds homology; unmatched ones are
    appended — these are typically the small tandem DUPs with no PE signal.
    """
    out = list(primary)
    for junc in recovered:
        matched = False
        for existing in out:
            if (
                existing.svtype == junc.svtype
                and abs(existing.start - junc.start) <= tol
                and abs(existing.end - junc.end) <= tol
            ):
                existing.start, existing.end = junc.start, junc.end
                existing.cipos = existing.ciend = (0, 0)
                existing.sr_start = max(existing.sr_start, junc.sr_start)
                existing.sr_end = max(existing.sr_end, junc.sr_end)
                existing.cig = max(existing.cig, junc.cig)
                existing.precise = True
                if junc.homlen >= 0:
                    existing.homlen = junc.homlen
                matched = True
                break
        if not matched:
            out.append(junc)
    out.sort(key=lambda j: (j.start, j.end))
    return out
