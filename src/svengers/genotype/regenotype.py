"""Re-genotype a sample's own calls from its evidence store.

The BAF enters a call twice at calling time: in the genotype likelihood
(GT, GQ) and as the filter features BAFN / BAFDEV. Once ``svengers baf`` has
rewritten a store's BAF (no QUAL threshold), the calls still carry the old
values. ``regenotype_calls`` recomputes GT, GQ, CN, CNQ, RD, BAFN and BAFDEV
of every call of the raw per-sample VCF from the store, with the call's own
junction support (PE, SR, CG) and source kept as they are, so the filter can
be re-run and a model trained on features that match the store.

INSERT-source calls are genotyped by their detector, not by the likelihood,
and are copied unchanged.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pysam

from ..cohort.genotyper import StoreGenotyper
from ..config import Config
from ..extract.store import EvidenceStore
from ..variants import SOURCE_INSERT, CNVCall
from ..vcfout import span, spec_coords

logger = logging.getLogger(__name__)

_FIELDS = ("GT", "GQ", "CN", "CNQ", "RD", "BAFN", "BAFDEV")


def regenotype_calls(vcf_in: str | Path, store_path: str | Path, vcf_out: str | Path,
                     config: Config | None = None) -> dict[str, int]:
    """Write ``vcf_out``: ``vcf_in`` with the store-derived fields recomputed.
    Returns counts: calls, regenotyped, unchanged (INSERT / no depth data),
    gt_changed."""
    store = EvidenceStore(store_path)
    sg = StoreGenotyper(store, config or Config(), rewritten_baf=True)
    counts = {"calls": 0, "regenotyped": 0, "unchanged": 0, "gt_changed": 0}
    with pysam.VariantFile(str(vcf_in)) as vin:
        spec = spec_coords(vin.header)
        hdr = vin.header.copy()
        hdr.add_line(f"##svengers_regenotyped=baf_source:{store.meta.get('baf_source') or 'unknown'}")
        with pysam.VariantFile(str(vcf_out), "wz", header=hdr) as vout:
            for rec in vin:
                counts["calls"] += 1
                rec.translate(hdr)
                svtype = rec.info.get("SVTYPE")
                source = rec.info.get("SOURCE") or ""
                data = sg._contig_data(rec.contig) if svtype in ("DEL", "DUP") else None
                if source == SOURCE_INSERT or data is None or data.ploidy == 0:
                    counts["unchanged"] += 1
                    vout.write(rec)
                    continue
                start, end = span(rec, spec)
                fmt = rec.samples[0]
                call = CNVCall(contig=rec.contig, start=start, end=end, svtype=svtype,
                               source=source, pe_support=int(fmt.get("PE") or 0),
                               sr_support=int(fmt.get("SR") or 0),
                               cig_support=int(fmt.get("CG") or 0))
                sg.genotyper.genotype(call, data)
                old_gt = fmt.get("GT")
                new_gt = write_genotype(fmt, call)
                counts["regenotyped"] += 1
                counts["gt_changed"] += tuple(old_gt or ()) != tuple(new_gt or (None, None))
                vout.write(rec)
    pysam.tabix_index(str(vcf_out), preset="vcf", force=True)
    return counts


def write_genotype(fmt, call: CNVCall):
    """Write the store-derived fields of a genotyped call (``_FIELDS``) into a
    record's sample; returns the new GT (None when undetermined)."""
    new_gt = None if call.gt == (-1, -1) else call.gt
    fmt["GT"] = new_gt if new_gt is not None else (None, None)
    fmt["GQ"] = call.gq
    if call.cn >= 0:
        fmt["CN"] = call.cn
        fmt["CNQ"] = round(call.cn_posterior, 4)
    if not np.isnan(call.mean_ratio):
        fmt["RD"] = round(call.mean_ratio, 4)
    fmt["BAFN"] = call.baf_n
    if np.isnan(call.baf_dev):
        if "BAFDEV" in fmt and fmt.get("BAFDEV") is not None:
            del fmt["BAFDEV"]
    else:
        fmt["BAFDEV"] = round(call.baf_dev, 4)
    return new_gt
