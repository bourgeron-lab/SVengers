"""Genomic-context annotation of calls: segdup/LCR/exclude overlap, GC, popCovStdev."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..reference import read_bed

logger = logging.getLogger(__name__)


class IntervalSet:
    """Merged sorted intervals per contig with fast overlap-fraction queries."""

    def __init__(self, bed: dict[str, np.ndarray]):
        self.starts: dict[str, np.ndarray] = {}
        self.ends: dict[str, np.ndarray] = {}
        for contig, ivs in bed.items():
            if len(ivs) == 0:
                continue
            ivs = ivs[np.argsort(ivs[:, 0])]
            starts, ends = [], []
            cur_s, cur_e = int(ivs[0, 0]), int(ivs[0, 1])
            for s, e in ivs[1:]:
                if s <= cur_e:
                    cur_e = max(cur_e, int(e))
                else:
                    starts.append(cur_s)
                    ends.append(cur_e)
                    cur_s, cur_e = int(s), int(e)
            starts.append(cur_s)
            ends.append(cur_e)
            self.starts[contig] = np.asarray(starts, dtype=np.int64)
            self.ends[contig] = np.asarray(ends, dtype=np.int64)

    @classmethod
    def from_bed(cls, path: str | Path) -> "IntervalSet":
        return cls(read_bed(path))

    def overlap_fraction(self, contig: str, start: int, end: int) -> float:
        if end <= start or contig not in self.starts:
            return 0.0
        starts, ends = self.starts[contig], self.ends[contig]
        i = int(np.searchsorted(ends, start, side="right"))
        total = 0
        while i < len(starts) and starts[i] < end:
            total += min(end, ends[i]) - max(start, starts[i])
            i += 1
        return total / (end - start)


@dataclass
class ContextAnnotator:
    segdup: IntervalSet | None = None
    lcr: IntervalSet | None = None
    exclude: IntervalSet | None = None
    popcov_bw: object | None = None    # opened pyBigWig handle

    @classmethod
    def from_paths(
        cls,
        segdup_bed: str = "",
        lcr_bed: str = "",
        exclude_bed: str = "",
        pop_cov_stdev_bw: str = "",
    ) -> "ContextAnnotator":
        ann = cls()
        if segdup_bed:
            ann.segdup = IntervalSet.from_bed(segdup_bed)
        if lcr_bed:
            ann.lcr = IntervalSet.from_bed(lcr_bed)
        if exclude_bed:
            ann.exclude = IntervalSet.from_bed(exclude_bed)
        if pop_cov_stdev_bw:
            try:
                import pyBigWig

                ann.popcov_bw = pyBigWig.open(pop_cov_stdev_bw)
            except Exception:
                logger.warning("Could not open popCovStdev bigWig %s", pop_cov_stdev_bw)
        return ann

    def annotate(self, contig: str, start: int, end: int) -> dict[str, float]:
        out: dict[str, float] = {}
        out["segdup_ovl"] = (
            self.segdup.overlap_fraction(contig, start, end) if self.segdup else 0.0
        )
        out["lcr_ovl"] = self.lcr.overlap_fraction(contig, start, end) if self.lcr else 0.0
        out["exclude_ovl"] = (
            self.exclude.overlap_fraction(contig, start, end) if self.exclude else 0.0
        )
        if self.popcov_bw is not None:
            try:
                val = self.popcov_bw.stats(contig, start, min(end, self.popcov_bw.chroms(contig)))
                out["popcov_stdev"] = float(val[0]) if val and val[0] is not None else 0.0
            except Exception:
                out["popcov_stdev"] = 0.0
        else:
            out["popcov_stdev"] = 0.0
        return out
