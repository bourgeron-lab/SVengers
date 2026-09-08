"""Trio-based specificity metrics as a function of the quality threshold.

For each threshold, genotypes below it (on GQ, or SQ when the cohort VCF
carries it) are treated as no-calls, then the Mendelian-violation rate and
the het transmission rate are recomputed. MVR falling and transmission
approaching 50% as the threshold rises is the headline truth-set-free
specificity figure.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pysam

from ..cohort.mendel import _gt_of, classify_trio
from ..ped import Pedigree

logger = logging.getLogger(__name__)

SWEEP_COLUMNS = [
    "threshold", "field", "n_genotyped", "n_carrier",
    "n_violation", "mvr", "n_transmit_opportunity", "transmission", "n_denovo",
]


def _qual_of(rec, sample: str, field: str) -> float:
    fmt = rec.samples.get(sample)
    if fmt is None:
        return 0.0
    val = fmt.get(field)
    return float(val) if val is not None else 0.0


def sweep(
    cohort_vcf: str | Path,
    pedigree: Pedigree,
    thresholds: list[float] = (0, 3, 7, 10, 13, 17, 20, 30),
    field: str = "GQ",
    out_tsv: str | Path | None = None,
) -> list[dict]:
    trios = pedigree.trios()
    with pysam.VariantFile(str(cohort_vcf)) as vcf:
        vcf_samples = set(vcf.header.samples)
        usable = [t for t in trios if all(i.sample_id in vcf_samples for i in t)]
        if not usable:
            raise ValueError("No complete trio present in the cohort VCF")
        if field not in vcf.header.formats:
            raise ValueError(f"FORMAT/{field} absent from {cohort_vcf}")
        # load once: per site, per trio, (gts, quals)
        site_data = []
        for rec in vcf:
            per_trio = []
            for child, father, mother in usable:
                trio_gt = []
                for ind in (child, father, mother):
                    gt, _gq = _gt_of(rec, ind.sample_id)
                    trio_gt.append((gt, _qual_of(rec, ind.sample_id, field)))
                per_trio.append(trio_gt)
            site_data.append(per_trio)

    rows = []
    for th in thresholds:
        n_geno = n_carrier = n_viol = n_opp = n_trans = n_denovo = 0
        for per_trio in site_data:
            for (cgt, cq), (fgt, fq), (mgt, mq) in per_trio:
                if cgt is None or fgt is None or mgt is None:
                    continue
                if min(cq, fq, mq) < th:
                    continue
                n_geno += 1
                inheritance, violation = classify_trio(cgt, fgt, mgt)
                if sum(cgt) > 0:
                    n_carrier += 1
                if violation:
                    n_viol += 1
                if inheritance == "denovo":
                    n_denovo += 1
                for pgt in (fgt, mgt):
                    if sum(pgt) == 1:
                        n_opp += 1
                        if sum(cgt) > 0:
                            n_trans += 1
        rows.append({
            "threshold": th, "field": field,
            "n_genotyped": n_geno, "n_carrier": n_carrier,
            "n_violation": n_viol,
            "mvr": round(n_viol / n_geno, 5) if n_geno else "",
            "n_transmit_opportunity": n_opp,
            "transmission": round(n_trans / n_opp, 4) if n_opp else "",
            "n_denovo": n_denovo,
        })
        logger.info(
            "%s>=%s: %d genotyped, MVR %s, transmission %s, %d de novo",
            field, th, n_geno, rows[-1]["mvr"], rows[-1]["transmission"], n_denovo,
        )

    if out_tsv:
        with open(out_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=SWEEP_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    return rows
