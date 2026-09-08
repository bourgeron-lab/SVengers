"""VCF 4.3 output (symbolic <DEL>/<DUP> records) via pysam."""

from __future__ import annotations

import math
from pathlib import Path

import pysam

from . import __version__
from .reference import Contig
from .variants import CNVCall

FILTER_DEFS = [
    ("MinSize", "Call smaller than the minimum reported size"),
    ("MaxN", "Call overlaps excluded/N regions above threshold"),
    ("LowMap", "Call lies in low-mappability/segdup context without split-read support"),
    ("NoEvidence", "Insufficient supporting evidence"),
    ("PopStdev", "Depth-only call in a population high-variance region"),
    ("LowSQ", "Calibrated quality below threshold"),
    ("CohortArtifact", "Recurrent cohort artifact locus"),
    ("LowPosterior", "Depth-only call with CN posterior below threshold"),
]

INFO_DEFS = [
    ("SVTYPE", "1", "String", "Type of structural variant"),
    ("END", "1", "Integer", "End position of the variant"),
    ("SVLEN", "1", "Integer", "Signed length of the variant"),
    ("CIPOS", "2", "Integer", "Confidence interval around POS"),
    ("CIEND", "2", "Integer", "Confidence interval around END"),
    ("SOURCE", "1", "String", "Evidence source: DEPTH, JUNCTION or BOTH"),
    ("NBINS", "1", "Integer", "Number of depth bins spanned"),
    ("MASKFRAC", "1", "Float", "Fraction of spanned bins that are masked"),
    ("IMPRECISE", "0", "Flag", "Breakpoints not base-pair resolved"),
    ("HOMLEN", "1", "Integer", "Reference micro-homology length at the junction"),
]

FORMAT_DEFS = [
    ("GT", "1", "String", "Genotype"),
    ("GQ", "1", "Integer", "Genotype quality (phred)"),
    ("SQ", "1", "Float", "Calibrated variant quality from the Vedette filter (phred)"),
    ("CN", "1", "Integer", "Most likely integer copy number"),
    ("CNQ", "1", "Float", "Copy-number posterior probability"),
    ("RD", "1", "Float", "Mean GC-corrected copy ratio over the call"),
    ("PE", "1", "Integer", "Supporting discordant pairs"),
    ("SR", "1", "Integer", "Supporting split/clipped reads"),
    ("BAFN", "1", "Integer", "Informative het sites within the call"),
    ("BAFDEV", "1", "Float", "Mean absolute deviation of het VAF from 0.5"),
]


def make_header(contigs: list[Contig], sample: str, command: str = "") -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_line(f"##source=vedette-{__version__}")
    if command:
        header.add_line(f"##vedette_command={command}")
    for c in contigs:
        header.contigs.add(c.name, length=c.length)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    for name, desc in FILTER_DEFS:
        header.filters.add(name, None, None, desc)
    for name, number, vtype, desc in INFO_DEFS:
        header.info.add(name, number, vtype, desc)
    for name, number, vtype, desc in FORMAT_DEFS:
        header.formats.add(name, number, vtype, desc)
    header.add_sample(sample)
    return header


def _fmt_gt(gt: tuple[int, int]) -> tuple:
    return (None, None) if gt == (-1, -1) else gt


def write_vcf(
    calls: list[CNVCall],
    contigs: list[Contig],
    sample: str,
    path: str | Path,
    command: str = "",
) -> None:
    header = make_header(contigs, sample, command)
    order = {c.name: i for i, c in enumerate(contigs)}
    calls = sorted(calls, key=lambda c: (order.get(c.contig, 1_000_000), c.start))

    path = str(path)
    mode = "wz" if path.endswith(".gz") else "w"
    with pysam.VariantFile(path, mode, header=header) as out:
        for i, call in enumerate(calls):
            rec = out.new_record(
                contig=call.contig,
                start=call.start,
                stop=call.end,
                alleles=("N", f"<{call.svtype}>"),
                id=f"vedette_{call.svtype}_{i + 1}",
            )
            rec.info["SVTYPE"] = call.svtype
            rec.info["SVLEN"] = call.svlen
            rec.info["SOURCE"] = call.source
            rec.info["NBINS"] = call.n_bins
            rec.info["MASKFRAC"] = round(call.masked_frac, 4)
            if not call.precise:
                rec.info["IMPRECISE"] = True
                rec.info["CIPOS"] = call.cipos
                rec.info["CIEND"] = call.ciend
            for key, value in call.info_extra.items():
                rec.info[key] = value
            for f in call.filters or ["PASS"]:
                rec.filter.add(f)

            fmt = rec.samples[sample]
            fmt["GT"] = _fmt_gt(call.gt)
            fmt["GQ"] = int(call.gq)
            if not math.isnan(call.sq):
                fmt["SQ"] = round(float(call.sq), 2)
            if call.cn >= 0:
                fmt["CN"] = call.cn
                fmt["CNQ"] = round(call.cn_posterior, 4)
            if not math.isnan(call.mean_ratio):
                fmt["RD"] = round(call.mean_ratio, 4)
            fmt["PE"] = call.pe_support
            fmt["SR"] = call.sr_support
            fmt["BAFN"] = call.baf_n
            if not math.isnan(call.baf_dev):
                fmt["BAFDEV"] = round(call.baf_dev, 4)
            out.write(rec)
    if mode == "wz":
        pysam.tabix_index(path, preset="vcf", force=True)
