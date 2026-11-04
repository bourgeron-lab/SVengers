"""Truth-set labeling of Vedette calls via Truvari.

Runs `truvari bench` (unfiltered permissive calls vs a GIAB-style truth VCF,
restricted to the confident-region BED), then joins the TP/FP status back onto
the per-call feature matrix produced by `vedette filter --features-out`.

Calls outside the confident regions get no label and are dropped from
training. DUP calls are matched against insertion-represented truth via
Truvari's --dup-to-ins.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

from ..vcfout import span, spec_coords

logger = logging.getLogger(__name__)

TRUVARI_ARGS = [
    "--pctseq", "0",        # symbolic calls carry no sequence
    "--pctsize", "0.7",
    "--refdist", "1000",
    "--sizemin", "50",
    "--sizemax", "20000000",
    "--dup-to-ins",
]


def _truvari_bin() -> str:
    exe = Path(sys.executable).with_name("truvari")
    return str(exe) if exe.exists() else "truvari"


def run_truvari_bench(
    calls_vcf: str | Path,
    truth_vcf: str | Path,
    include_bed: str | Path,
    reference: str | Path | None,
    out_dir: str | Path,
    extra_args: list[str] | None = None,
) -> Path:
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)  # truvari refuses an existing directory
    out_dir.parent.mkdir(parents=True, exist_ok=True)  # truvari only mkdirs the leaf
    cmd = [
        _truvari_bin(), "bench",
        "-b", str(truth_vcf),
        "-c", str(calls_vcf),
        "--includebed", str(include_bed),
        *(["--reference", str(reference)] if reference else []),
        "-o", str(out_dir),
        *TRUVARI_ARGS,
        *(extra_args or []),
    ]
    logger.info("truvari: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "truvari bench failed:\n"
            f"--- stderr ---\n{result.stderr[-2000:]}\n"
            f"--- stdout ---\n{result.stdout[-2000:]}"
        )
    return out_dir


def _ids_of(vcf_path: Path) -> set[str]:
    ids: set[str] = set()
    if not vcf_path.exists():
        return ids
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            if rec.id:
                ids.add(rec.id)
    return ids


def _fp_intervals(vcf_path: str | Path):
    """IntervalSet of FP-labeled calls (for parent-mode recurrent-artifact labels)."""
    import numpy as np

    from ..filter.annotate import IntervalSet

    bed: dict[str, list[tuple[int, int]]] = {}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        spec = spec_coords(vcf.header)
        for rec in vcf:
            bed.setdefault(rec.contig, []).append(span(rec, spec))
    return IntervalSet(
        {c: np.array(sorted(ivs), dtype=np.int64) for c, ivs in bed.items()}
    )


def label_features(
    features_tsv: str | Path,
    truvari_dir: str | Path,
    out_tsv: str | Path,
    parent_mode: bool = False,
    recurrent_fp_vcf: str | Path | None = None,
    min_fp_overlap: float = 0.5,
) -> dict[str, int]:
    """Join truvari TP/FP status onto the feature matrix -> labeled TSV.

    ``parent_mode`` implements transmission-based semi-supervision for samples
    WITHOUT their own truth set (HG003/HG004 vs the HG002 truth): a truth match
    is a certain TP (transmitted variant); a non-match is labeled FP only when
    it overlaps a known recurrent-FP locus (``recurrent_fp_vcf``, e.g. the
    proband's FP calls) — otherwise it may be a private true variant and is
    dropped as ambiguous.
    """
    truvari_dir = Path(truvari_dir)
    tp_ids = _ids_of(truvari_dir / "tp-comp.vcf.gz")
    fp_ids = _ids_of(truvari_dir / "fp.vcf.gz")
    overlap = tp_ids & fp_ids
    if overlap:
        raise RuntimeError(f"{len(overlap)} call IDs are both TP and FP — duplicate IDs?")
    fp_loci = _fp_intervals(recurrent_fp_vcf) if parent_mode and recurrent_fp_vcf else None
    if parent_mode and fp_loci is None:
        raise ValueError("parent_mode requires recurrent_fp_vcf")

    counts = {"tp": 0, "fp": 0, "ambiguous": 0, "unlabeled": 0}
    with open(features_tsv, newline="") as fin, open(out_tsv, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        writer.writerow([*header, "label"])
        for row in reader:
            call_id = row[0]
            if call_id in tp_ids:
                writer.writerow([*row, 1])
                counts["tp"] += 1
            elif call_id in fp_ids:
                if not parent_mode:
                    writer.writerow([*row, 0])
                    counts["fp"] += 1
                else:
                    contig = row[idx["contig"]]
                    start, end = int(row[idx["start"]]), int(row[idx["end"]])
                    if fp_loci.overlap_fraction(contig, start, end) >= min_fp_overlap:
                        writer.writerow([*row, 0])
                        counts["fp"] += 1
                    else:
                        counts["ambiguous"] += 1  # possibly a private true variant
            else:
                counts["unlabeled"] += 1  # outside confident regions / size range
    logger.info(
        "labels: %(tp)d TP, %(fp)d FP, %(ambiguous)d ambiguous (dropped), "
        "%(unlabeled)d unlabeled (dropped)", counts,
    )
    return counts
