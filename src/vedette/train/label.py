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
        raise RuntimeError(f"truvari bench failed:\n{result.stderr[-3000:]}")
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


def label_features(
    features_tsv: str | Path,
    truvari_dir: str | Path,
    out_tsv: str | Path,
) -> dict[str, int]:
    """Join truvari TP/FP status onto the feature matrix -> labeled TSV."""
    truvari_dir = Path(truvari_dir)
    tp_ids = _ids_of(truvari_dir / "tp-comp.vcf.gz")
    fp_ids = _ids_of(truvari_dir / "fp.vcf.gz")
    overlap = tp_ids & fp_ids
    if overlap:
        raise RuntimeError(f"{len(overlap)} call IDs are both TP and FP — duplicate IDs?")

    counts = {"tp": 0, "fp": 0, "unlabeled": 0}
    with open(features_tsv, newline="") as fin, open(out_tsv, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        writer.writerow([*header, "label"])
        for row in reader:
            call_id = row[0]
            if call_id in tp_ids:
                writer.writerow([*row, 1])
                counts["tp"] += 1
            elif call_id in fp_ids:
                writer.writerow([*row, 0])
                counts["fp"] += 1
            else:
                counts["unlabeled"] += 1  # outside confident regions / size range
    logger.info(
        "labels: %(tp)d TP, %(fp)d FP, %(unlabeled)d unlabeled (dropped)", counts
    )
    return counts
