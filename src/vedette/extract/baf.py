"""B-allele fraction evidence.

Fast path: a tabix-indexed per-site VAF bedgraph (GHFC ``{ID}.vaf.bedgraph.gz``:
chrom, pos, pos+1, vaf — no depth column, dp is stored as 0 = unknown).
VCF path: a single-sample SNV VCF (DeepVariant, GATK) read with
``read_vaf_vcf``, and ``backfill_baf`` to rewrite the BAF of existing stores
from it.

The GHFC track keeps ``PASS && QUAL >= 30 && DP >= 10``. DeepVariant 1.9
gives many true heterozygotes a QUAL below 30: on chr20 an ordinary 39x
genome keeps 44k of its 56k PASS heterozygotes (DP >= 10), and a 75x old
(IP) genome 1.6k of 32k, the ones it loses having a median VAF of 0.50. No
QUAL threshold here.
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
    path: str, contig: str, min_dp: int = 8, pass_only: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(pos, vaf, dp) from a single-sample SNV VCF with AD/DP FORMAT fields:
    biallelic SNVs with a depth (DP, else the sum of AD) of at least
    ``min_dp``, PASS (or
    unfiltered ".") with ``pass_only``. DeepVariant writes its reference calls
    as RefCall records, at a VAF of 0.1-0.3 for most of them: they are no
    heterozygous sites."""
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
            if pass_only and any(f != "PASS" for f in rec.filter.keys()):
                continue
            sample = rec.samples[0]
            ad = sample.get("AD")
            if ad is None or ad[0] is None:
                continue
            total = sum(a for a in ad if a is not None)
            # the site depth is FORMAT/DP when given (the GHFC track's
            # DP >= 10), else the allele depths' sum
            depth = sample.get("DP") if "DP" in sample else None
            depth = total if depth is None else depth
            if depth < min_dp or not total:
                continue
            pos.append(rec.pos - 1)  # store 0-based
            vaf.append(ad[1] / total)
            dp.append(min(depth, 65535))
    return (
        np.asarray(pos, np.uint32),
        np.asarray(vaf, np.float32),
        np.asarray(dp, np.uint16),
    )


def backfill_baf(store_path, vcf_path: str, min_dp: int = 10,
                 force: bool = False) -> tuple[int, int, bool]:
    """(sites, heterozygous sites, written): rewrite the BAF of an existing
    evidence store from a single-sample SNV VCF (``read_vaf_vcf``, PASS,
    depth >= ``min_dp``), without the CRAM. The previous arrays are
    kept as ``<contig>/baf_track/`` (once: a second backfill does not
    overwrite them); the store moves to schema 0.6 with ``baf_source``.

    Only the store changes: the BAF features of the sample's calls (BAFN,
    BAFDEV) were computed at calling time and stay as they are."""
    import shutil
    import time
    from pathlib import Path

    import zarr

    from .store import SCHEMA_VERSION, EvidenceStore, write_baf_group

    store_path = Path(store_path)
    store = EvidenceStore(store_path)
    source = f"vcf:{Path(vcf_path).name} PASS DP>={min_dp}"
    if (store.meta.get("baf_source") or "").startswith("vcf:") and not force:
        return 0, 0, False
    # read everything first: a VCF with no site for the store's contigs
    # (another naming, another sample) must not wipe its BAF
    new = {c: read_vaf_vcf(str(vcf_path), c, min_dp) for c in store.contigs}
    n_sites = sum(len(p) for p, _, _ in new.values())
    if not n_sites:
        raise ValueError(f"{vcf_path}: no PASS SNV with DP >= {min_dp} on the store's contigs")
    root = zarr.open_group(str(store_path), mode="r+")
    n_het = 0
    for contig, (pos, vaf, dp) in new.items():
        cdir = store_path / contig
        if (cdir / "baf").is_dir():
            if (cdir / "baf_track").exists():
                shutil.rmtree(cdir / "baf")
            else:
                (cdir / "baf").rename(cdir / "baf_track")
        if len(pos):
            write_baf_group(root[contig], pos, vaf, dp)
        n_het += int(((vaf > 0.1) & (vaf < 0.9)).sum())
    root.attrs.update(schema_version=SCHEMA_VERSION, baf_source=source,
                      baf_backfilled=time.strftime("%Y-%m-%dT%H:%M:%S"))
    return n_sites, n_het, True
