"""Pedigree annotation of a cohort VCF: inheritance, Mendelian violations,
de novo candidates, transmission rate.

The Mendelian-violation rate and the transmission rate (expected ~50% for a
true heterozygous variant) are *free specificity metrics*: they quantify the
false-positive/genotyping-error load without any truth set.

chrX and chrY follow the sex of the child: a son's chrX comes from his
mother (his genotype is hemizygous, 0/0 or 1/1), a daughter gets her
father's only chrX, chrY goes from father to son; in the pseudo-autosomal
regions only "carrier child, both parents 0/0" counts. The child's sex is
the ploidy its own genotypes show on chrX (reference CN 1 or 2), the
pedigree sex as a fallback. Read as autosomes, every son 1/1 on chrX with a
0/1 mother and a 0/0 father was a violation.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..ped import Pedigree
from ..vcfout import span, spec_coords

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


# GRCh38 pseudo-autosomal regions, 0-based half-open
PAR = {"chrX": ((10_000, 2_781_479), (155_701_382, 156_030_895)),
       "chrY": ((10_000, 2_781_479), (56_887_902, 57_217_415))}


def site_mode(contig: str, start: int, stop: int) -> str:
    """'auto', 'x', 'y' or 'par' (a site overlapping a pseudo-autosomal region)."""
    c = contig if contig.startswith("chr") else f"chr{contig}"
    if c not in PAR:
        return "auto"
    if any(start < e and stop > b for b, e in PAR[c]):
        return "par"
    return "x" if c == "chrX" else "y"


def transmit_parents(mode: str, father_gt, mother_gt) -> tuple:
    """Parent genotypes that give a ~50 % transmission opportunity: none from
    a hemizygous father (his chrX always goes to his daughters, his chrY to
    his sons)."""
    if mode == "x":
        return (mother_gt,)
    if mode == "y":
        return ()
    return (father_gt, mother_gt)


def classify_trio(
    child_gt, father_gt, mother_gt, mode: str = "auto", child_male: bool = False,
) -> tuple[str, bool]:
    """(inheritance label, is_violation) for a biallelic site; ``mode`` from
    site_mode."""
    c = sum(child_gt)
    f = sum(father_gt)
    m = sum(mother_gt)
    if mode == "x" and child_male:
        if c == 0:
            return "none", m == 2
        return ("denovo", True) if m == 0 else ("maternal", False)
    if mode == "y":
        if not child_male or c == 0:
            return "none", child_male and f > 0
        return ("denovo", True) if f == 0 else ("paternal", False)
    if mode == "x":      # daughter: the father's only chrX
        fa = f > 0
        violation = (
            (c == 2 and (not fa or m == 0))
            or (c == 1 and ((not fa and m == 0) or (fa and m == 2)))
            or (c == 0 and (fa or m == 2))
        )
    elif mode == "par":
        violation = c > 0 and f == 0 and m == 0
    else:
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


def hemizygous_children(
    cohort_vcf: str | Path, trios, pedigree: Pedigree | None = None,
    max_records: int = 3000,
) -> set[str]:
    """Children with one chrX: median CN of their 0/0 genotypes on chrX
    outside the PARs (every few records, up to ``max_records``) below 1.5;
    the pedigree sex when the VCF has no chrX copy number."""
    import numpy as np

    children = sorted({t[0].sample_id for t in trios})
    cns: dict[str, list[int]] = {s: [] for s in children}
    with pysam.VariantFile(str(cohort_vcf)) as vcf:
        contig = next((c for c in ("chrX", "X") if c in vcf.header.contigs), None)
        if contig and "CN" in vcf.header.formats and children:
            vcf.subset_samples(children)
            try:
                recs = list(vcf.fetch(contig))
            except ValueError:
                recs = []
            spec = spec_coords(vcf.header)
            recs = [r for r in recs if site_mode(contig, *span(r, spec)) == "x"]
            for rec in recs[:: max(1, len(recs) // max_records)]:
                for s in children:
                    fmt = rec.samples[s]
                    gt, cn = fmt.get("GT"), fmt.get("CN")
                    if gt is not None and cn is not None and all(a == 0 for a in gt):
                        cns[s].append(int(cn))
    male = set()
    for s in children:
        if len(cns[s]) >= 10:
            if float(np.median(cns[s])) < 1.5:
                male.add(s)
        elif pedigree is not None and pedigree.sex_of(s) == 1:
            male.add(s)
    return male


def _ped_region(task) -> tuple[dict[str, list[int]], str]:
    """One region of the cohort VCF: per-child counters (n_genotyped,
    n_child_carrier, n_denovo, n_violation, n_transmit_opportunity,
    n_transmitted) and the report rows written to ``part`` (no header)."""
    from ..parallel import fetch_region

    (cohort_vcf, region, usable, part,
     min_denovo_child_gq, min_denovo_parent_gq, male) = task
    counters = {child.sample_id: [0, 0, 0, 0, 0, 0] for child, _f, _m in usable}
    with pysam.VariantFile(str(cohort_vcf)) as vcf, open(part, "w", newline="") as fh:
        # decode only the trio members' genotypes
        vcf.subset_samples(sorted({ind.sample_id for t in usable for ind in t}))
        writer = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t")
        records = fetch_region(vcf, region) if region is not None else vcf
        spec = spec_coords(vcf.header)
        for rec in records:
            start, end = span(rec, spec)
            mode = site_mode(rec.contig, start, end)
            for child, father, mother in usable:
                cgt, cgq = _gt_of(rec, child.sample_id)
                fgt, fgq = _gt_of(rec, father.sample_id)
                mgt, mgq = _gt_of(rec, mother.sample_id)
                if cgt is None or fgt is None or mgt is None:
                    continue
                st = counters[child.sample_id]
                st[0] += 1
                inheritance, violation = classify_trio(
                    cgt, fgt, mgt, mode, child.sample_id in male
                )
                if sum(cgt) > 0:
                    st[1] += 1
                if violation:
                    st[3] += 1
                # transmission: each het parent gives one opportunity
                for pgt in transmit_parents(mode, fgt, mgt):
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
                    # 0-based half-open, as a BED interval
                    "start": start, "end": end,
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
    male = hemizygous_children(cohort_vcf, usable, pedigree)
    logger.info("ped: %d children with one chrX", len(male))
    regions = split_regions(cohort_vcf, jobs)
    tmp = tempfile.mkdtemp(prefix=".ped_parts_", dir=Path(out_report).parent)
    try:
        tasks = [
            (str(cohort_vcf), region, usable, f"{tmp}/part{k:05d}.tsv",
             min_denovo_child_gq, min_denovo_parent_gq, male)
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
