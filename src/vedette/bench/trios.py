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
    pass_only: bool = False,
    exclude_rescue: bool = False,
) -> list[dict]:
    trios = pedigree.trios()
    thresholds = sorted(thresholds)
    # single streaming pass: per (site, trio) update every threshold <= min GQ
    # (storing all trio genotypes first costs tens of GB at 2676 x 80k sites)
    counters = [
        {"n_geno": 0, "n_carrier": 0, "n_viol": 0, "n_opp": 0, "n_trans": 0, "n_denovo": 0}
        for _ in thresholds
    ]
    with pysam.VariantFile(str(cohort_vcf)) as vcf:
        vcf_samples = set(vcf.header.samples)
        usable = [t for t in trios if all(i.sample_id in vcf_samples for i in t)]
        if not usable:
            raise ValueError("No complete trio present in the cohort VCF")
        if field not in vcf.header.formats:
            raise ValueError(f"FORMAT/{field} absent from {cohort_vcf}")
        for rec in vcf:
            if pass_only and set(rec.filter.keys()) - {"PASS"}:
                continue
            if exclude_rescue and rec.info.get("RESCUE"):
                continue
            for child, father, mother in usable:
                cgt, _ = _gt_of(rec, child.sample_id)
                fgt, _ = _gt_of(rec, father.sample_id)
                mgt, _ = _gt_of(rec, mother.sample_id)
                if cgt is None or fgt is None or mgt is None:
                    continue
                minq = min(
                    _qual_of(rec, ind.sample_id, field)
                    for ind in (child, father, mother)
                )
                inheritance, violation = classify_trio(cgt, fgt, mgt)
                carrier = sum(cgt) > 0
                opp = sum(1 for pgt in (fgt, mgt) if sum(pgt) == 1)
                for k, th in enumerate(thresholds):
                    if th > minq:
                        break
                    c = counters[k]
                    c["n_geno"] += 1
                    if carrier:
                        c["n_carrier"] += 1
                        c["n_trans"] += opp
                    if violation:
                        c["n_viol"] += 1
                    if inheritance == "denovo":
                        c["n_denovo"] += 1
                    c["n_opp"] += opp

    rows = []
    for th, c in zip(thresholds, counters, strict=True):
        n_geno, n_carrier = c["n_geno"], c["n_carrier"]
        n_viol, n_opp, n_trans, n_denovo = (
            c["n_viol"], c["n_opp"], c["n_trans"], c["n_denovo"]
        )
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
