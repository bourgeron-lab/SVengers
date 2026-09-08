"""Reference genome handling: contigs, binning, GC/N tracks, exclusion masks.

Per-reference companion tracks (GC fraction and N fraction per bin) are computed
once and cached to an ``.npz`` in a user-writable resources directory, keyed by
bin size; every sample then reuses them.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------- contigs


@dataclass(frozen=True)
class Contig:
    name: str
    length: int
    kind: str  # "autosome" | "X" | "Y" | "MT" | "other"


def _contig_kind(name: str) -> str:
    base = name.removeprefix("chr")
    if base.isdigit():
        return "autosome"
    if base == "X":
        return "X"
    if base == "Y":
        return "Y"
    if base in ("M", "MT"):
        return "MT"
    return "other"


def read_fai(fasta: str | Path) -> list[Contig]:
    """Parse the .fai next to ``fasta`` into an ordered contig list."""
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(f"FASTA index not found: {fai} (run `samtools faidx`)")
    contigs = []
    with open(fai) as fh:
        for line in fh:
            fields = line.split("\t")
            contigs.append(Contig(fields[0], int(fields[1]), _contig_kind(fields[0])))
    return contigs


class GenomeBins:
    """Fixed-size binning of a set of contigs."""

    def __init__(self, contigs: list[Contig], bin_size: int):
        self.bin_size = bin_size
        self.contigs = {c.name: c for c in contigs}
        self.n_bins = {c.name: -(-c.length // bin_size) for c in contigs}

    def bin_slice(self, contig: str, start: int, end: int) -> slice:
        """Bins fully or partially overlapping [start, end) (0-based, half-open)."""
        end = min(end, self.contigs[contig].length)
        return slice(start // self.bin_size, -(-end // self.bin_size))

    def bin_starts(self, contig: str) -> np.ndarray:
        return np.arange(self.n_bins[contig], dtype=np.int64) * self.bin_size


# ------------------------------------------------------------------------- GC/N tracks

_GC_BYTES = frozenset(b"GCgc")
_N_BYTES = frozenset(b"Nn")


def _seq_tracks(seq: bytes, bin_size: int) -> tuple[np.ndarray, np.ndarray]:
    """GC fraction (of non-N bases) and N fraction per bin for one contig sequence."""
    arr = np.frombuffer(seq, dtype=np.uint8)
    n_bins = -(-len(arr) // bin_size)
    padded = np.zeros(n_bins * bin_size, dtype=np.uint8)
    padded[: len(arr)] = arr
    padded = padded.reshape(n_bins, bin_size)

    is_gc = np.isin(padded, np.frombuffer(b"GCgc", dtype=np.uint8))
    is_n = np.isin(padded, np.frombuffer(b"Nn", dtype=np.uint8))
    # padding bytes (0) count as N so trailing partial bins are masked
    is_n |= padded == 0

    n_count = is_n.sum(axis=1)
    denom = np.maximum(bin_size - n_count, 1)
    gc = is_gc.sum(axis=1) / denom
    return gc.astype(np.float32), (n_count / bin_size).astype(np.float32)


def reference_id(fasta: str | Path) -> str:
    """Cheap stable identifier of a reference: sha1 of its .fai content."""
    with open(str(fasta) + ".fai", "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()[:16]


def tracks_path(
    resources_dir: str | Path, fasta: str | Path, bin_size: int, contigs: list[Contig]
) -> Path:
    contig_key = hashlib.sha1("|".join(c.name for c in contigs).encode()).hexdigest()[:8]
    return Path(resources_dir) / (
        f"reftracks_{reference_id(fasta)}_{bin_size}_{contig_key}.npz"
    )


def build_tracks(
    fasta: str | Path,
    contigs: list[Contig],
    bin_size: int,
    resources_dir: str | Path,
) -> Path:
    """Compute and cache GC/N per-bin tracks for the given contigs."""
    import pysam

    out = tracks_path(resources_dir, fasta, bin_size, contigs)
    if out.exists():
        return out
    Path(resources_dir).mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    with pysam.FastaFile(str(fasta)) as fa:
        for contig in contigs:
            logger.info("GC/N track: %s", contig.name)
            seq = fa.fetch(contig.name).encode()
            gc, nfrac = _seq_tracks(seq, bin_size)
            arrays[f"gc/{contig.name}"] = gc
            arrays[f"nfrac/{contig.name}"] = nfrac
    tmp = out.with_suffix(".npz.tmp.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.rename(out)
    logger.info("Cached reference tracks: %s", out)
    return out


class ReferenceTracks:
    """Lazy reader of the cached GC/N tracks."""

    def __init__(self, path: str | Path):
        self._npz = np.load(path)

    def gc(self, contig: str) -> np.ndarray:
        return self._npz[f"gc/{contig}"]

    def nfrac(self, contig: str) -> np.ndarray:
        return self._npz[f"nfrac/{contig}"]


# ---------------------------------------------------------------------------- BED masks


def read_bed(path: str | Path) -> dict[str, np.ndarray]:
    """Read a (possibly gzipped) BED into {contig: array[(start, end)]} (sorted)."""
    intervals: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return {
        contig: np.array(sorted(ivs), dtype=np.int64) for contig, ivs in intervals.items()
    }


def bed_bin_overlap(
    bed: dict[str, np.ndarray], contig: str, n_bins: int, bin_size: int
) -> np.ndarray:
    """Fraction of each bin covered by BED intervals (float32, len n_bins)."""
    cover = np.zeros(n_bins, dtype=np.float32)
    ivs = bed.get(contig)
    if ivs is None or len(ivs) == 0:
        return cover
    starts = np.clip(ivs[:, 0], 0, n_bins * bin_size)
    ends = np.clip(ivs[:, 1], 0, n_bins * bin_size)
    keep = ends > starts
    starts, ends = starts[keep], ends[keep]
    first = starts // bin_size
    last = (ends - 1) // bin_size  # inclusive

    single = first == last
    np.add.at(cover, first[single], (ends[single] - starts[single]) / bin_size)

    multi = ~single
    np.add.at(cover, first[multi], ((first[multi] + 1) * bin_size - starts[multi]) / bin_size)
    np.add.at(cover, last[multi], (ends[multi] - last[multi] * bin_size) / bin_size)
    # full middle bins (intervals spanning >2 bins are rare; loop only over those)
    for s, e in zip(first[multi], last[multi], strict=True):
        if e > s + 1:
            cover[s + 1 : e] += 1.0
    return np.minimum(cover, 1.0)


def build_bin_mask(
    n_bins: int,
    bin_size: int,
    contig: str,
    nfrac: np.ndarray,
    exclude_bed: dict[str, np.ndarray] | None,
    max_nfrac: float = 0.1,
    max_exclude_frac: float = 0.5,
) -> np.ndarray:
    """Boolean mask of bins to EXCLUDE from depth modeling."""
    mask = nfrac > max_nfrac
    if exclude_bed is not None:
        mask |= bed_bin_overlap(exclude_bed, contig, n_bins, bin_size) > max_exclude_frac
    return mask
