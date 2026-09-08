"""Pedigree annotation of a cohort VCF: inheritance, Mendelian violations,
de novo candidates, transmission rate.

The Mendelian-violation rate and the transmission rate (expected ~50% for a
true heterozygous variant) are *free specificity metrics*: they quantify the
false-positive/genotyping-error load without any truth set.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from ..ped import Pedigree

logger = logging.getLogger(__name__)

REPORT_COLUMNS = [
    "site_id", "contig", "start", "end", "svtype", "family_id", "child",
    "child_gt", "child_gq", "father_gt", "father_gq", "mother_gt", "mother_gq",
    "inheritance", "violation", "dnq",
]


@dataclass
class TrioStats:
    family_id: str
    child: str
    n_genotyped: int = 0
    n_child_carrier: int = 0
    n_denovo: int = 0
    n_violation: int = 0
    n_transmit_opportunity: int = 0   # parent het-carrier, child genotyped
    n_transmitted: int = 0

    @property
    def violation_rate(self) -> float:
        return self.n_violation / self.n_genotyped if self.n_genotyped else 0.0

    @property
    def transmission_rate(self) -> float:
        return (
            self.n_transmitted / self.n_transmit_opportunity
            if self.n_transmit_opportunity else float("nan")
        )


def _gt_of(rec, sample: str) -> tuple[tuple[int, ...] | None, int]:
    fmt = rec.samples.get(sample)
    if fmt is None:
        return None, 0
    gt = fmt.get("GT")
    if gt is None or any(a is None for a in gt):
        return None, 0
    return tuple(gt), int(fmt.get("GQ") or 0)


def classify_trio(
    child_gt, father_gt, mother_gt
) -> tuple[str, bool]:
    """(inheritance label, is_violation) for a biallelic site."""
    c = sum(child_gt)
    f = sum(father_gt)
    m = sum(mother_gt)
    violation = (
        (c == 2 and (f == 0 or m == 0))
        or (c == 0 and (f == 2 or m == 2))
        or (c == 1 and f == 0 and m == 0)   # het child, both parents hom-ref
        or (c == 0 and f == 2 and m == 2)
    )
    if c == 0:
        return "none", violation
    if f == 0 and m == 0:
        return "denovo", True
    if f > 0 and m > 0:
        return "biparental", violation
    if f > 0:
        return "paternal", violation
    return "maternal", violation


def run_ped_annotate(
    cohort_vcf: str | Path,
    pedigree: Pedigree,
    out_report: str | Path,
    out_summary: str | Path | None = None,
    min_denovo_child_gq: int = 30,
    min_denovo_parent_gq: int = 20,
) -> list[TrioStats]:
    trios = pedigree.trios()
    stats = {
        t[0].sample_id: TrioStats(t[0].family_id, t[0].sample_id) for t in trios
    }
    with (
        pysam.VariantFile(str(cohort_vcf)) as vcf,
        open(out_report, "w", newline="") as fh,
    ):
        vcf_samples = set(vcf.header.samples)
        usable = [
            t for t in trios
            if all(ind.sample_id in vcf_samples for ind in t)
        ]
        logger.info("ped: %d/%d trios present in the VCF", len(usable), len(trios))
        writer = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for rec in vcf:
            for child, father, mother in usable:
                cgt, cgq = _gt_of(rec, child.sample_id)
                fgt, fgq = _gt_of(rec, father.sample_id)
                mgt, mgq = _gt_of(rec, mother.sample_id)
                if cgt is None or fgt is None or mgt is None:
                    continue
                st = stats[child.sample_id]
                st.n_genotyped += 1
                inheritance, violation = classify_trio(cgt, fgt, mgt)
                if sum(cgt) > 0:
                    st.n_child_carrier += 1
                if violation:
                    st.n_violation += 1
                # transmission: each het parent gives one opportunity
                for pgt in (fgt, mgt):
                    if sum(pgt) == 1:
                        st.n_transmit_opportunity += 1
                        if sum(cgt) > 0:
                            st.n_transmitted += 1
                dnq = ""
                if inheritance == "denovo":
                    if cgq >= min_denovo_child_gq and min(fgq, mgq) >= min_denovo_parent_gq:
                        st.n_denovo += 1
                        dnq = str(min(cgq, fgq, mgq))
                    else:
                        inheritance = "denovo_lowconf"
                if sum(cgt) == 0 and not violation:
                    continue  # only report informative rows
                writer.writerow({
                    "site_id": rec.id, "contig": rec.contig,
                    "start": rec.start, "end": rec.stop,
                    "svtype": rec.info.get("SVTYPE", ""),
                    "family_id": child.family_id, "child": child.sample_id,
                    "child_gt": "/".join(map(str, cgt)), "child_gq": cgq,
                    "father_gt": "/".join(map(str, fgt)), "father_gq": fgq,
                    "mother_gt": "/".join(map(str, mgt)), "mother_gq": mgq,
                    "inheritance": inheritance,
                    "violation": int(violation), "dnq": dnq,
                })

    result = [stats[t[0].sample_id] for t in usable]
    if out_summary:
        with open(out_summary, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "family_id", "child", "n_genotyped", "n_child_carrier",
                "n_denovo", "n_violation", "violation_rate",
                "transmission_rate",
            ])
            for st in result:
                writer.writerow([
                    st.family_id, st.child, st.n_genotyped, st.n_child_carrier,
                    st.n_denovo, st.n_violation,
                    f"{st.violation_rate:.4f}",
                    f"{st.transmission_rate:.4f}" if st.n_transmit_opportunity else "NA",
                ])
    for st in result:
        logger.info(
            "trio %s/%s: %d genotyped, %d carrier, %d de novo, MVR %.3f, transmission %.3f",
            st.family_id, st.child, st.n_genotyped, st.n_child_carrier,
            st.n_denovo, st.violation_rate,
            st.transmission_rate if st.n_transmit_opportunity else float("nan"),
        )
    return result
