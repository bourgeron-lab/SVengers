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

De novo and the child's own call. The cohort genotyper force-genotypes every
site in every sample, so a child can come out a carrier at a site other
samples nominated while the caller found nothing in that child. A true de
novo is private to the child: the child should hold a call there. Given the
per-sample filtered VCFs (``calls_dir``), a de novo whose child has no call of
the same type at the site (reciprocal overlap >= ``own_overlap``, any FILTER)
is labelled ``denovo_genotype_only`` and does not count in ``n_denovo``; the
report's ``child_call`` column gives the FILTER of the child's call (``none``
when absent). On dataset_paper (665 trios), 90 % of the PASS GQ >= 20 de novo
had no call of the child; in a blind image review of 40 of them, none of the
26 without a call was a de novo at the site (3 were fragments of a larger
event the child's own call covers), while 10 of 14 with a PASS call of the
child, no parental call and a private site were.

A call of the child is not enough. On dataset_paper again, 90 % of the de
novo that kept a call of the child at a PASS site held a LowSQ call only, and
most sat on sites that came out de novo in several unrelated children: a
genotyping miss in the parents at a polymorphic or artifact-prone site, not a
new mutation. With ``calls_dir`` a de novo is therefore confirmed
(``denovo``, counted in ``n_denovo``) only when the three properties the image
review found in the true ones hold:

* the child's own call at the site is PASS;
* neither parent holds a call of the same type there (any FILTER: a LowSQ
  parental call is the evidence the forced 0/0 missed);
* few samples outside the child's family carry the site in the cohort
  genotypes: at most ``max_other_carriers`` or ``max_other_carrier_frac`` of
  the pedigree's samples, whichever is larger (10 on dataset_paper's 2 000).
  A new mutation is private to its family, but a recurrent NAHR event is
  not: two unrelated children of dataset_paper carry the same de novo 15q11-13
  BP1-BP3 duplication, and the forced genotypes of a large event add a few
  carriers; with at most 2 other carriers, every de novo 15q11-13 and 22q11
  duplication of dataset_paper was lost.

Other GQ-confident de novo with a call of the child are ``denovo_weak``
(``n_denovo_weak``). The report gives the parents' own calls
(``parent_call``, father,mother) and the carriers outside the family
(``other_carriers``). On dataset_paper this takes the confirmed DUP from 680
(v12b) and 2 563 (v14b) to ~40, the same events for both models: mostly large
private events (an autosomal trisomy shows up as one call per fragment).
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
    "inheritance", "violation", "dnq", "child_call", "parent_call", "other_carriers",
]


@dataclass
class TrioStats:
    family_id: str
    child: str
    n_genotyped: int = 0
    n_child_carrier: int = 0
    n_denovo: int = 0
    n_violation: int = 0
    n_denovo_genotype_only: int = 0   # high-confidence de novo without a call of the child
    n_transmit_opportunity: int = 0   # parent het-carrier, child genotyped
    n_transmitted: int = 0
    n_denovo_weak: int = 0            # child call not PASS, parental call or shared site

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


class OwnCalls:
    """The child's own calls at a site, from its filtered single-sample VCF.

    A region worker walks the sites in order and checks every child at each,
    so a child's carrier calls are loaded once per contig (one tabix fetch)
    and the cache is dropped when the contig changes."""

    def __init__(self, paths: dict[str, str], min_overlap: float = 0.5):
        self.paths = paths
        self.min_overlap = min_overlap
        self._contig: str | None = None
        self._cache: dict[str, dict] = {}

    def _load(self, sample: str, contig: str) -> dict:
        import numpy as np

        if contig != self._contig:
            self._contig, self._cache = contig, {}
        got = self._cache.get(sample)
        if got is not None:
            return got
        by: dict[str, list] = {}
        with pysam.VariantFile(self.paths[sample]) as vf:
            spec, name = spec_coords(vf.header), list(vf.header.samples)[0]
            try:
                recs = vf.fetch(contig)
            except ValueError:
                recs = []
            for rec in recs:
                t = rec.info.get("SVTYPE")
                gt = rec.samples[name].get("GT")
                if t is None or gt is None or not any(a == 1 for a in gt if a is not None):
                    continue
                s, e = span(rec, spec)
                flt = ";".join(rec.filter.keys()) or "PASS"
                by.setdefault(t, []).append((s, e, flt, float(rec.samples[name].get("SQ") or 0.0)))
        got = {t: (np.array([x[0] for x in v]), np.array([x[1] for x in v]),
                   [x[2] for x in v], np.array([x[3] for x in v])) for t, v in by.items()}
        self._cache[sample] = got
        return got

    def status(self, sample: str, contig: str, start: int, end: int, svtype: str) -> str | None:
        """FILTER of the child's best matching carrier call ('PASS' first, then
        the highest SQ), 'none' without one, None when the sample has no VCF."""
        import numpy as np

        if sample not in self.paths:
            return None
        x = self._load(sample, contig).get(svtype)
        if x is None:
            return "none"
        s, e, flt, sq = x
        ov = np.minimum(e, end) - np.maximum(s, start)
        m = (ov > 0) & (ov >= self.min_overlap * (end - start)) & (ov >= self.min_overlap * (e - s))
        idx = np.flatnonzero(m)
        if not len(idx):
            return "none"
        return flt[min(idx, key=lambda i: (flt[i] != "PASS", -sq[i]))]

    def close(self) -> None:
        self._contig, self._cache = None, {}


def own_call_paths(calls_dir: str | Path, samples) -> dict[str, str]:
    """{sample: <calls_dir>/<S>/<S>.filtered.vcf.gz} for those present (Nextflow layout)."""
    out = {}
    for s in samples:
        p = Path(calls_dir) / s / f"{s}.filtered.vcf.gz"
        if p.exists():
            out[s] = str(p)
    return out


def _carriers_by_family(rec, family_of: dict[str, str]) -> dict[str, int]:
    """Carriers (a non-reference allele) of the record per pedigree family,
    over the decoded samples."""
    out: dict[str, int] = {}
    for s, fmt in rec.samples.items():
        gt = fmt.get("GT")
        if gt and any(a for a in gt if a is not None):
            f = family_of.get(s)
            if f is not None:
                out[f] = out.get(f, 0) + 1
    return out


def _ped_region(task) -> tuple[dict[str, list[int]], str]:
    """One region of the cohort VCF: per-child counters (n_genotyped,
    n_child_carrier, n_denovo, n_violation, n_transmit_opportunity,
    n_transmitted, n_denovo_genotype_only, n_denovo_weak) and the report rows
    written to ``part`` (no header)."""
    from ..parallel import fetch_region

    (cohort_vcf, region, usable, part, min_denovo_child_gq, min_denovo_parent_gq,
     male, own_paths, own_overlap, family_of, max_other_carriers) = task
    own = OwnCalls(own_paths, own_overlap) if own_paths else None
    counters = {child.sample_id: [0, 0, 0, 0, 0, 0, 0, 0] for child, _f, _m in usable}
    with pysam.VariantFile(str(cohort_vcf)) as vcf, open(part, "w", newline="") as fh:
        # decode only the samples of the pedigree's families
        vcf.subset_samples(sorted(s for s in vcf.header.samples if s in family_of))
        writer = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t")
        records = fetch_region(vcf, region) if region is not None else vcf
        spec = spec_coords(vcf.header)
        for rec in records:
            start, end = span(rec, spec)
            mode = site_mode(rec.contig, start, end)
            carriers_by_family = None   # {family: n carriers}, on the record's first de novo
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
                dnq = child_call = parent_call = other = ""
                if inheritance == "denovo":
                    confident = cgq >= min_denovo_child_gq and min(fgq, mgq) >= min_denovo_parent_gq
                    svtype = rec.info.get("SVTYPE", "")
                    if own is not None:
                        child_call = own.status(
                            child.sample_id, rec.contig, start, end, svtype
                        ) or ""
                    if own is not None and confident and child_call not in ("none", ""):
                        pc = [own.status(p.sample_id, rec.contig, start, end, svtype) or "NA"
                              for p in (father, mother)]
                        parent_call = ",".join(pc)
                        if carriers_by_family is None:
                            carriers_by_family = _carriers_by_family(rec, family_of)
                        n_other = sum(carriers_by_family.values()) - carriers_by_family.get(
                            child.family_id, 0)
                        other = str(n_other)
                        strict = (
                            child_call.split(";")[0] == "PASS"
                            and all(c in ("none", "NA") for c in pc)
                            and n_other <= max_other_carriers
                        )
                    else:
                        strict = True
                    if not confident:
                        inheritance = "denovo_lowconf"
                    elif child_call == "none":
                        inheritance = "denovo_genotype_only"
                        st[6] += 1
                    elif not strict:
                        inheritance = "denovo_weak"
                        st[7] += 1
                    else:
                        st[2] += 1
                        dnq = str(min(cgq, fgq, mgq))
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
                    "violation": int(violation), "dnq": dnq, "child_call": child_call,
                    "parent_call": parent_call, "other_carriers": other,
                })
    if own is not None:
        own.close()
    return counters, part


def run_ped_annotate(
    cohort_vcf: str | Path,
    pedigree: Pedigree,
    out_report: str | Path,
    out_summary: str | Path | None = None,
    min_denovo_child_gq: int = 30,
    min_denovo_parent_gq: int = 20,
    jobs: int = 1,
    calls_dir: str | Path | None = None,
    own_overlap: float = 0.5,
    max_other_carriers: int = 2,
    max_other_carrier_frac: float = 0.005,
) -> list[TrioStats]:
    """``jobs`` > 1 (indexed VCF) annotates position-ordered regions in
    parallel and joins their rows in file order: same report as serial.
    ``calls_dir`` (``<dir>/<S>/<S>.filtered.vcf.gz``) checks each de novo
    against the child's own calls (module docstring)."""
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
    own_paths: dict[str, str] = {}
    if calls_dir is not None:
        members = sorted({ind.sample_id for t in usable for ind in t})
        own_paths = own_call_paths(calls_dir, members)
        logger.info("ped: own calls of %d/%d trio members under %s", len(own_paths), len(members), calls_dir)
    family_of = {s: ind.family_id for s, ind in pedigree.individuals.items()}
    n_ped = sum(1 for s in vcf_samples if s in family_of)
    max_other = max(max_other_carriers, int(max_other_carrier_frac * n_ped))
    if calls_dir is not None:
        logger.info("ped: a confirmed de novo site has <= %d carriers outside the family", max_other)
    regions = split_regions(cohort_vcf, jobs)
    tmp = tempfile.mkdtemp(prefix=".ped_parts_", dir=Path(out_report).parent)
    try:
        tasks = [
            (str(cohort_vcf), region, usable, f"{tmp}/part{k:05d}.tsv",
             min_denovo_child_gq, min_denovo_parent_gq, male, own_paths, own_overlap,
             family_of, max_other)
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
                st.n_denovo_genotype_only += c[6]
                st.n_denovo_weak += c[7]
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
                "transmission_rate", "n_denovo_genotype_only", "n_denovo_weak",
            ])
            for st in result:
                writer.writerow([
                    st.family_id, st.child, st.n_genotyped, st.n_child_carrier,
                    st.n_denovo, st.n_violation,
                    f"{st.violation_rate:.4f}",
                    f"{st.transmission_rate:.4f}" if st.n_transmit_opportunity else "NA",
                    st.n_denovo_genotype_only if own_paths else "NA",
                    st.n_denovo_weak if own_paths else "NA",
                ])
    for st in result:
        logger.info(
            "trio %s/%s: %d genotyped, %d carrier, %d de novo, MVR %.3f, transmission %.3f",
            st.family_id, st.child, st.n_genotyped, st.n_child_carrier,
            st.n_denovo, st.violation_rate,
            st.transmission_rate if st.n_transmit_opportunity else float("nan"),
        )
    return result
