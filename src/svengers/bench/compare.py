"""Multi-caller comparison on the same truth set.

Harmonizes external caller outputs (symbolic SV VCFs like ClinSV/Manta/DELLY/
dysgu, or WisecondorX aberration BEDs) into normalized DEL/DUP VCFs, then runs
the same stratified Truvari evaluation on each — including SVengers itself
before and after filtering, to isolate the filter's marginal value.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pysam

from ..reference import Contig
from ..vcfout import add_coords_line, span, spec_coords, sv_record
from .truth import Stratifier, evaluate

logger = logging.getLogger(__name__)


def _norm_header(contigs: list[Contig]) -> pysam.VariantHeader:
    h = pysam.VariantHeader()
    add_coords_line(h)
    for c in contigs:
        h.contigs.add(c.name, length=c.length)
    h.add_line('##ALT=<ID=DEL,Description="Deletion">')
    h.add_line('##ALT=<ID=DUP,Description="Duplication">')
    h.info.add("SVTYPE", "1", "String", "SV type")
    h.info.add("SVLEN", "1", "Integer", "Signed SV length")
    h.info.add("END", "1", "Integer", "End position")
    h.formats.add("GT", "1", "String", "Genotype")
    h.add_sample("S")
    return h


def _write_norm(
    out_vcf: Path, contigs: list[Contig], records: list[tuple[str, int, int, str]], caller: str
) -> int:
    """records: (contig, start, end, svtype), 0-based half-open."""
    order = {c.name: i for i, c in enumerate(contigs)}
    known = set(order)
    records = sorted(
        (r for r in records if r[0] in known and r[2] > r[1]),
        key=lambda r: (order[r[0]], r[1]),
    )
    header = _norm_header(contigs)
    with pysam.VariantFile(str(out_vcf), "wz", header=header) as out:
        for i, (contig, start, end, svtype) in enumerate(records):
            rec = sv_record(out, contig, start, end, svtype, f"{caller}_{i + 1}")
            rec.samples["S"]["GT"] = (0, 1)
            out.write(rec)
    pysam.tabix_index(str(out_vcf), preset="vcf", force=True)
    return len(records)


def normalize_sv_vcf(
    in_vcf: str | Path, out_vcf: str | Path, contigs: list[Contig],
    caller: str, pass_only: bool = True, min_size: int = 50,
) -> int:
    """Extract DEL/DUP from any symbolic SV VCF into the normalized form."""
    records = []
    with pysam.VariantFile(str(in_vcf)) as vcf:
        spec = spec_coords(vcf.header)
        for rec in vcf:
            svtype = str(rec.info.get("SVTYPE", ""))
            if svtype.startswith("DUP"):
                svtype = "DUP"
            if svtype not in ("DEL", "DUP"):
                continue
            if pass_only and set(rec.filter.keys()) - {"PASS", "."}:
                continue
            start, end = span(rec, spec)
            if end - start < min_size:
                continue
            records.append((rec.contig, start, end, svtype))
    n = _write_norm(Path(out_vcf), contigs, records, caller)
    logger.info("%s: %d DEL/DUP calls normalized from %s", caller, n, in_vcf)
    return n


def wisecondorx_to_vcf(
    aberrations_bed: str | Path, out_vcf: str | Path, contigs: list[Contig],
    caller: str = "wisecondorx",
) -> int:
    """Convert a WisecondorX *_aberrations(.chr).bed (loss/gain) to the normalized VCF."""
    records = []
    with open(aberrations_bed) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        gain_loss_col = idx.get("source", len(header) - 1)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            contig = fields[idx.get("chr", 0)]
            if not contig.startswith("chr"):
                contig = f"chr{contig}"
            kind = fields[gain_loss_col].strip().lower()
            if kind not in ("gain", "loss"):
                continue
            records.append((
                contig,
                int(fields[idx.get("start", 1)]),
                int(fields[idx.get("end", 2)]),
                "DUP" if kind == "gain" else "DEL",
            ))
    n = _write_norm(Path(out_vcf), contigs, records, caller)
    logger.info("%s: %d aberrations converted from %s", caller, n, aberrations_bed)
    return n


def run_compare(
    entries: list[tuple[str, str, str]],   # (name, kind: vcf|vcf-all|wisecondorx, path)
    truth_vcf: str | Path,
    include_bed: str | Path,
    contigs: list[Contig],
    out_dir: str | Path,
    reference: str | Path | None = None,
    stratifier: Stratifier | None = None,
) -> list[dict]:
    """Evaluate every caller with the same protocol; combined tidy table."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for name, kind, path in entries:
        norm = out_dir / f"{name}.norm.vcf.gz"
        if kind == "wisecondorx":
            n = wisecondorx_to_vcf(path, norm, contigs, caller=name)
        else:
            # "vcf-all" keeps non-PASS records (pre-thresholded operating points)
            n = normalize_sv_vcf(
                path, norm, contigs, caller=name, pass_only=(kind != "vcf-all")
            )
        if n == 0:
            logger.warning("%s: no calls after normalization, skipping", name)
            continue
        rows = evaluate(
            norm, truth_vcf, include_bed, out_dir / name,
            reference=reference, stratifier=stratifier,
        )
        for row in rows:
            all_rows.append({"caller": name, **row})

    out_tsv = out_dir / "comparison.tsv"
    if all_rows:
        with open(out_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(all_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)
        logger.info("comparison table -> %s", out_tsv)
    return all_rows
