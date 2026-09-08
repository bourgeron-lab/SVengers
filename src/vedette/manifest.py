"""Cohort sample manifest: sample_id -> CRAM path (+ optional batch, sex).

TSV with header: sample_id, cram, [batch], [sex], [depth_bedgraph], [vaf_bedgraph].
``discover`` auto-builds one from a known directory layout (currently the GHFC
layout: samples/{last_char}/{second_to_last_char}/{ID}/sequences/{ID}.<ref>.cram
with optional precomputed depth/VAF bedgraphs alongside).
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

_COLUMNS = ("sample_id", "cram", "batch", "sex", "depth_bedgraph", "vaf_bedgraph")


@dataclass(frozen=True)
class Sample:
    sample_id: str
    cram: str
    batch: str = ""
    sex: int = 0
    depth_bedgraph: str = ""
    vaf_bedgraph: str = ""


@dataclass
class Manifest:
    samples: dict[str, Sample] = field(default_factory=dict)

    def add(self, sample: Sample) -> None:
        if sample.sample_id in self.samples:
            raise ValueError(f"Duplicate sample_id {sample.sample_id!r}")
        self.samples[sample.sample_id] = sample

    def write(self, path: str | Path) -> None:
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(_COLUMNS)
            for s in self.samples.values():
                writer.writerow(
                    [s.sample_id, s.cram, s.batch, s.sex, s.depth_bedgraph, s.vaf_bedgraph]
                )

    @classmethod
    def read(cls, path: str | Path) -> Manifest:
        manifest = cls()
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"sample_id", "cram"}
            if not required <= set(reader.fieldnames or ()):
                raise ValueError(f"{path}: manifest needs at least columns {sorted(required)}")
            for row in reader:
                manifest.add(
                    Sample(
                        sample_id=row["sample_id"].strip(),
                        cram=row["cram"].strip(),
                        batch=(row.get("batch") or "").strip(),
                        sex=int(row.get("sex") or 0),
                        depth_bedgraph=(row.get("depth_bedgraph") or "").strip(),
                        vaf_bedgraph=(row.get("vaf_bedgraph") or "").strip(),
                    )
                )
        return manifest

    def validate(self) -> list[str]:
        """Return a list of problems (missing CRAM/index files)."""
        problems = []
        for s in self.samples.values():
            cram = Path(s.cram)
            if not cram.exists():
                problems.append(f"{s.sample_id}: CRAM not found: {cram}")
            elif not (
                cram.with_suffix(cram.suffix + ".crai").exists()
                or cram.with_suffix(".crai").exists()
            ):
                problems.append(f"{s.sample_id}: no .crai index for {cram}")
        return problems


def discover_ghfc(
    root: str | Path,
    sample_ids: list[str] | None = None,
    cram_suffix: str = ".GRCh38_GIABv3.cram",
) -> Manifest:
    """Build a manifest from the GHFC samples/{c1}/{c2}/{ID}/sequences/ layout.

    If ``sample_ids`` is given (e.g. from a pedigree), only those are looked up
    and missing ones are reported; otherwise the whole tree is walked.
    """
    root = Path(root)
    manifest = Manifest()

    def _add_from_dir(sample_id: str) -> bool:
        seqdir = root / sample_id[-1] / sample_id[-2] / sample_id / "sequences"
        cram = seqdir / f"{sample_id}{cram_suffix}"
        if not cram.exists():
            return False
        depth_bg = seqdir / f"{sample_id}.by1000.bedgraph.gz"
        vaf_bg = seqdir / f"{sample_id}.vaf.bedgraph.gz"
        manifest.add(
            Sample(
                sample_id=sample_id,
                cram=str(cram),
                depth_bedgraph=str(depth_bg) if depth_bg.exists() else "",
                vaf_bedgraph=str(vaf_bg) if vaf_bg.exists() else "",
            )
        )
        return True

    if sample_ids is not None:
        missing = [sid for sid in sample_ids if not _add_from_dir(sid)]
        for sid in missing:
            logger.warning("discover: no CRAM found for %s", sid)
    else:
        for shard1 in sorted(p for p in root.iterdir() if p.is_dir()):
            for shard2 in sorted(p for p in shard1.iterdir() if p.is_dir()):
                for sample_dir in sorted(p for p in shard2.iterdir() if p.is_dir()):
                    _add_from_dir(sample_dir.name)
    logger.info("discover: %d samples with CRAMs under %s", len(manifest.samples), root)
    return manifest
