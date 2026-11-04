"""VCF 4.3 output (symbolic <DEL>/<DUP> records) via pysam, and the one
conversion between VCF and internal coordinates.

Internally an event is 0-based half-open ``[start, end)``: the affected
bases. In the VCF (4.3 spec for symbolic alleles), POS is the padding base
before the event and END its last base, so POS = ``start`` and END = ``end``
in 1-based terms (pysam: ``rec.start = start - 1``, ``rec.stop = end``) and
SVLEN = +/-(end - start). Files carry ``##vedette_coords=padded``; Vedette
files without it (before v1.4.0) wrote POS = start + 1 and END = end + 1
(htslib sets END = POS + |SVLEN|), and are read with that convention.
"""

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
    ("CG", "1", "Integer", "Reads carrying the event as a CIGAR deletion"),
    ("IZ", "1", "Float", "Insert-size-shift z-score (INSERT-source calls)"),
    ("IP", "1", "Integer", "Spanning proper pairs behind the insert-shift estimate"),
    ("BAFN", "1", "Integer", "Informative het sites within the call"),
    ("BAFDEV", "1", "Float", "Mean absolute deviation of het VAF from 0.5"),
]


def make_header(contigs: list[Contig], sample: str, command: str = "") -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_line(f"##source=vedette-{__version__}")
    add_coords_line(header)
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


COORDS_KEY = "vedette_coords"
COORDS_VALUE = "padded"


def add_coords_line(header: pysam.VariantHeader) -> None:
    """Mark a header as written with the VCF-spec coordinates (module doc)."""
    if not any(r.key == COORDS_KEY for r in header.records):
        header.add_line(f"##{COORDS_KEY}={COORDS_VALUE}")


def spec_coords(header: pysam.VariantHeader) -> bool:
    """False only for a Vedette file written before v1.4.0 (``##source=vedette``
    without the coordinates line); any other VCF (truth sets, other callers)
    follows the spec."""
    records = list(header.records)
    if any(r.key == COORDS_KEY for r in records):
        return True
    return not any(r.key == "source" and str(r.value).startswith("vedette") for r in records)


def span(rec: pysam.VariantRecord, spec: bool = True) -> tuple[int, int]:
    """Internal 0-based half-open ``[start, end)`` of an SV record.

    Spec: POS is the padding base (``rec.start + 1``) and END the last base
    (``rec.stop``); this also holds for sequence-resolved deletions (REF =
    padding + deleted bases). Legacy Vedette files: ``[rec.start, rec.stop - 1)``.
    """
    if spec:
        return rec.start + 1, max(rec.stop, rec.start + 2)
    return rec.start, max(rec.stop - 1, rec.start + 1)


def sv_record(out, contig: str, start: int, end: int, svtype: str, rec_id: str | None = None):
    """New symbolic record for the internal span ``[start, end)``, with SVTYPE
    and SVLEN set (``out``: a VariantFile or VariantHeader). An event at the
    very contig start has no padding base: it is written from base 1."""
    start = max(start, 1)
    end = max(end, start + 1)
    rec = out.new_record(
        contig=contig, start=start - 1, stop=end, alleles=("N", f"<{svtype}>"), id=rec_id,
    )
    rec.info["SVTYPE"] = svtype
    # htslib recomputes END = POS + |SVLEN| from SVLEN: equal to ``end`` here
    rec.info["SVLEN"] = (end - start) * (1 if svtype == "DUP" else -1)
    return rec


def clamp_end(start: int, end: int, length: int | None) -> int:
    """Half-open end bounded by the contig, so the written END stays on it.

    A depth segment's last 1 kb bin overhangs the contig end (dataset_paper:
    24 cohort sites up to 1 kb past chrX/chr12/...).
    """
    if length is not None:
        end = min(end, length)
    return max(end, start + 1)


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
    lengths = {c.name: c.length for c in contigs}
    calls = sorted(calls, key=lambda c: (order.get(c.contig, 1_000_000), c.start))

    path = str(path)
    mode = "wz" if path.endswith(".gz") else "w"
    with pysam.VariantFile(path, mode, header=header) as out:
        for i, call in enumerate(calls):
            end = clamp_end(call.start, call.end, lengths.get(call.contig))
            rec = sv_record(out, call.contig, call.start, end, call.svtype,
                            f"vedette_{call.svtype}_{i + 1}")
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
            fmt["CG"] = call.cig_support
            if call.ins_z > 0:
                fmt["IZ"] = round(call.ins_z, 2)
                fmt["IP"] = call.ins_pairs
            fmt["BAFN"] = call.baf_n
            if not math.isnan(call.baf_dev):
                fmt["BAFDEV"] = round(call.baf_dev, 4)
            out.write(rec)
    if mode == "wz":
        pysam.tabix_index(path, preset="vcf", force=True)
