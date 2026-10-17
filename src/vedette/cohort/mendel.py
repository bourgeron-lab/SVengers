"""Pedigree annotation of a cohort VCF: inheritance, Mendelian violations,
de novo candidates, transmission rate.

The Mendelian-violation rate and the transmission rate (expected ~50% for a
true heterozygous variant) are *free specificity metrics*: they quantify the
false-positive/genotyping-error load without any truth set.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
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


def _ped_region(task) -> tuple[dict[str, list[int]], str]:
    """One region of the cohort VCF: per-child counters (n_genotyped,
    n_child_carrier, n_denovo, n_violation, n_transmit_opportunity,
    n_transmitted) and the report rows written to ``part`` (no header)."""
    from ..parallel import fetch_region

    (cohort_vcf, region, usable, part,
     min_denovo_child_gq, min_denovo_parent_gq) = task
    counters = {child.sample_id: [0, 0, 0, 0, 0, 0] for child, _f, _m in usable}
    with pysam.VariantFile(str(cohort_vcf)) as vcf, open(part, "w", newline="") as fh:
        # decode only the trio members' genotypes
        vcf.subset_samples(sorted({ind.sample_id for t in usable for ind in t}))
        writer = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t")
        records = fetch_region(vcf, region) if region is not None else vcf
        for rec in records:
            for child, father, mother in usable:
                cgt, cgq = _gt_of(rec, child.sample_id)
                fgt, fgq = _gt_of(rec, father.sample_id)
                mgt, mgq = _gt_of(rec, mother.sample_id)
                if cgt is None or fgt is None or mgt is None:
                    continue
                st = counters[child.sample_id]
                st[0] += 1
                inheritance, violation = classify_trio(cgt, fgt, mgt)
                if sum(cgt) > 0:
                    st[1] += 1
                if violation:
                    st[3] += 1
                # transmission: each het parent gives one opportunity
                for pgt in (fgt, mgt):
                    if sum(pgt) == 1:
                        st[4] += 1
                        if sum(cgt) > 0:
                            st[5] += 1
                dnq = ""
                if inheritance == "denovo":
                    if cgq >= min_denovo_child_gq and min(fgq, mgq) >= min_denovo_parent_gq:
                        st[2] += 1
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
    return counters, part


def run_ped_annotate(
    cohort_vcf: str | Path,
    pedigree: Pedigree,
    out_report: str | Path,
    out_summary: str | Path | None = None,
    min_denovo_child_gq: int = 30,
    min_denovo_parent_gq: int = 20,
    jobs: int = 1,
) -> list[TrioStats]:
    """``jobs`` > 1 (indexed VCF) annotates position-ordered regions in
    parallel and joins their rows in file order: same report as serial."""
    import shutil
    import tempfile

    from ..parallel import concat_text, pmap, split_regions

    trios = pedigree.trios()
    stats = {
        t[0].sample_id: TrioStats(t[0].family_id, t[0].sample_id) for t in trios
    }
    with pysam.VariantFile(str(cohort_vcf)) as vcf:
        vcf_samples = set(vcf.header.samples)
    usable = [
        t for t in trios
        if all(ind.sample_id in vcf_samples for ind in t)
    ]
    logger.info("ped: %d/%d trios present in the VCF", len(usable), len(trios))
    regions = split_regions(cohort_vcf, jobs)
    tmp = tempfile.mkdtemp(prefix=".ped_parts_", dir=Path(out_report).parent)
    try:
        tasks = [
            (str(cohort_vcf), region, usable, f"{tmp}/part{k:05d}.tsv",
             min_denovo_child_gq, min_denovo_parent_gq)
            for k, region in enumerate(regions)
        ]
        parts = [f"{tmp}/header.tsv"]
        with open(parts[0], "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t").writeheader()
        for counters, part in pmap(_ped_region, tasks, jobs):
            parts.append(part)
            for child, c in counters.items():
                st = stats[child]
                st.n_genotyped += c[0]
                st.n_child_carrier += c[1]
                st.n_denovo += c[2]
                st.n_violation += c[3]
                st.n_transmit_opportunity += c[4]
                st.n_transmitted += c[5]
        concat_text(parts, out_report)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

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
