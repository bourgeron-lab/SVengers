"""Stratified truth-set evaluation (Truvari-based).

Runs `truvari bench` then post-stratifies the four outcome sets (tp-base, fn
= truth side; tp-comp, fp = call side) by size class and genomic context, so
precision/recall are reported per stratum — what reviewers expect since the
Q100/CMRG era, instead of genome-wide averages that hide segdup collapse.

Truth SVs are classified DEL vs INS from their explicit alleles; DUP calls
are evaluated against INS-represented truth via --dup-to-ins upstream.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import pysam

from ..filter.annotate import IntervalSet
from ..train.label import run_truvari_bench

logger = logging.getLogger(__name__)

SIZE_CLASSES = [
    (50, 300, "50-300bp"),
    (300, 1_000, "300bp-1kb"),
    (1_000, 10_000, "1-10kb"),
    (10_000, 100_000, "10-100kb"),
    (100_000, 1_000_000, "100kb-1Mb"),
    (1_000_000, float("inf"), ">1Mb"),
]

EVAL_COLUMNS = [
    "svtype", "stratum_kind", "stratum",
    "tp_base", "fn", "recall", "tp_comp", "fp", "precision", "f1",
]


def _size_class(size: int) -> str:
    for lo, hi, name in SIZE_CLASSES:
        if lo <= size < hi:
            return name
    return "<50bp"


def _records(vcf_path: Path):
    """Yield (svtype, size, contig, start, end) for one truvari output VCF."""
    if not vcf_path.exists():
        return
    with pysam.VariantFile(str(vcf_path)) as vcf:
        # pysam raises on an INFO key the header does not declare instead of
        # returning None; CMRG's explicit-allele VCF declares neither
        has_svtype = "SVTYPE" in vcf.header.info
        has_svlen = "SVLEN" in vcf.header.info
        for rec in vcf:
            svtype = rec.info.get("SVTYPE") if has_svtype else None
            svlen = rec.info.get("SVLEN") if has_svlen else None
            if isinstance(svlen, tuple):
                svlen = svlen[0]
            if svtype in ("DEL", "DUP", "INS"):
                size = abs(int(svlen)) if svlen is not None else rec.stop - rec.start
            else:  # explicit-allele truth record
                rlen, alen = len(rec.ref), len(rec.alts[0]) if rec.alts else 0
                size = abs(alen - rlen)
                svtype = "DEL" if rlen > alen else "INS"
            if svtype == "DUP":
                svtype = "INS"  # matched in insertion space; report jointly
            yield svtype, size, rec.contig, rec.start, max(rec.stop, rec.start + 1)


class Stratifier:
    def __init__(self, segdup_bed: str = "", lcr_bed: str = "", tr_bed: str = ""):
        self.sets: list[tuple[str, IntervalSet]] = []
        for name, path in (("segdup", segdup_bed), ("lcr", lcr_bed), ("tandem_repeat", tr_bed)):
            if path:
                self.sets.append((name, IntervalSet.from_bed(path)))

    def context(self, contig: str, start: int, end: int) -> str:
        for name, ivs in self.sets:  # priority order: segdup > lcr > TR
            if ivs.overlap_fraction(contig, start, end) > 0.5:
                return name
        return "clean"


def evaluate(
    calls_vcf: str | Path,
    truth_vcf: str | Path,
    include_bed: str | Path,
    out_dir: str | Path,
    reference: str | Path | None = None,
    stratifier: Stratifier | None = None,
    truvari_extra: list[str] | None = None,
) -> list[dict]:
    """Run truvari + stratified aggregation; returns tidy rows (also written as TSV)."""
    out_dir = Path(out_dir)
    tdir = run_truvari_bench(
        calls_vcf, truth_vcf, include_bed, reference, out_dir / "truvari",
        extra_args=truvari_extra,
    )
    stratifier = stratifier or Stratifier()

    # outcome -> (side, weight into which counters)
    sets = {
        "tp_base": tdir / "tp-base.vcf.gz",
        "fn": tdir / "fn.vcf.gz",
        "tp_comp": tdir / "tp-comp.vcf.gz",
        "fp": tdir / "fp.vcf.gz",
    }
    counts: dict[tuple[str, str, str], dict[str, int]] = {}

    def bump(svtype: str, kind: str, stratum: str, outcome: str) -> None:
        key = (svtype, kind, stratum)
        counts.setdefault(key, {k: 0 for k in sets})[outcome] += 1

    for outcome, path in sets.items():
        for svtype, size, contig, start, end in _records(path):
            bump(svtype, "size", _size_class(size), outcome)
            bump(svtype, "context", stratifier.context(contig, start, end), outcome)
            bump(svtype, "all", "all", outcome)

    rows = []
    order = {name: i for i, (_l, _h, name) in enumerate(SIZE_CLASSES)}
    for (svtype, kind, stratum), c in sorted(
        counts.items(), key=lambda t: (t[0][0], t[0][1], order.get(t[0][2], 99), t[0][2])
    ):
        recall = c["tp_base"] / (c["tp_base"] + c["fn"]) if c["tp_base"] + c["fn"] else None
        precision = c["tp_comp"] / (c["tp_comp"] + c["fp"]) if c["tp_comp"] + c["fp"] else None
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision and recall and (precision + recall) else None
        )
        rows.append({
            "svtype": svtype, "stratum_kind": kind, "stratum": stratum,
            "tp_base": c["tp_base"], "fn": c["fn"],
            "recall": round(recall, 4) if recall is not None else "",
            "tp_comp": c["tp_comp"], "fp": c["fp"],
            "precision": round(precision, 4) if precision is not None else "",
            "f1": round(f1, 4) if f1 is not None else "",
        })

    out_tsv = out_dir / "eval_stratified.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=EVAL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = json.loads((tdir / "summary.json").read_text())
    logger.info(
        "truvari overall: precision %.4f recall %.4f f1 %.4f (TP %d, FP %d, FN %d)",
        summary.get("precision") or 0, summary.get("recall") or 0, summary.get("f1") or 0,
        summary.get("TP-comp", 0), summary.get("FP", 0), summary.get("FN", 0),
    )
    return rows
