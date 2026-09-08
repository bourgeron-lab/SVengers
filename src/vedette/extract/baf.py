"""B-allele fraction evidence.

Fast path: a tabix-indexed per-site VAF bedgraph (GHFC ``{ID}.vaf.bedgraph.gz``:
chrom, pos, pos+1, vaf — no depth column, dp is stored as 0 = unknown).
Self-contained path: pileup of the CRAM at sites of a user-supplied SNV VCF
(implemented with the CRAM scan; see scan.py TODO v1.1).
"""

from __future__ import annotations

import logging

import numpy as np
import pysam

logger = logging.getLogger(__name__)


def read_vaf_bedgraph(
    path: str, contig: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(pos, vaf, dp) arrays for one contig from a per-site VAF bedgraph."""
    pos: list[int] = []
    vaf: list[float] = []
    with pysam.TabixFile(path) as tbx:
        if contig not in tbx.contigs:
            return (
                np.array([], np.uint32),
                np.array([], np.float32),
                np.array([], np.uint16),
            )
        for row in tbx.fetch(contig, parser=pysam.asTuple()):
            pos.append(int(row[1]))
            vaf.append(float(row[3]))
    dp = np.zeros(len(pos), dtype=np.uint16)  # unknown depth in fast path
    return np.asarray(pos, np.uint32), np.asarray(vaf, np.float32), dp


def read_vaf_vcf(
    path: str, contig: str, min_dp: int = 8
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(pos, vaf, dp) from a single-sample SNV VCF with AD/DP FORMAT fields."""
    pos: list[int] = []
    vaf: list[float] = []
    dp: list[int] = []
    with pysam.VariantFile(path) as vcf:
        try:
            it = vcf.fetch(contig)
        except ValueError:
            return (
                np.array([], np.uint32),
                np.array([], np.float32),
                np.array([], np.uint16),
            )
        for rec in it:
            if len(rec.alleles) != 2 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue  # biallelic SNVs only
            sample = rec.samples[0]
            ad = sample.get("AD")
            if ad is None or ad[0] is None:
                continue
            total = sum(a for a in ad if a is not None)
            if total < min_dp:
                continue
            pos.append(rec.pos - 1)  # store 0-based
            vaf.append(ad[1] / total)
            dp.append(min(total, 65535))
    return (
        np.asarray(pos, np.uint32),
        np.asarray(vaf, np.float32),
        np.asarray(dp, np.uint16),
    )
