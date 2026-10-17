"""Process-level parallelism for the cohort stages (merge, cohort-filter, ped,
report): a cohort VCF cut into position-ordered regions, and an executor.

Every cohort stage ran on one core (99 % CPU, ~14 h end to end on the
2670-sample cohort while the genotyping step used 200 jobs). The work is
independent per region, so each stage maps over regions and stitches the
results back IN FILE ORDER: the output is identical to the serial run.

A record belongs to the region holding its START, so a call spanning a
region boundary is read once (``fetch`` alone returns it in both regions).

Executor: ProcessPoolExecutor, not multiprocessing.Pool — an OOM-killed
worker hangs Pool forever but raises BrokenProcessPool here (see
extract/scan.py). Fork context, so that large read-only arrays built by the
parent before the pool starts are shared copy-on-write; the parent must not
hold an open htslib handle at that point (handles must not cross fork).
"""

from __future__ import annotations

import logging
import os
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import pysam

REGION_SIZE = 10_000_000

logger = logging.getLogger(__name__)


def default_jobs() -> int:
    """CPUs this process may use (the SLURM allocation on a compute node)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:  # not Linux
        return max(1, os.cpu_count() or 1)


def resolve_jobs(jobs: int) -> int:
    """0 or less = every allowed CPU."""
    return jobs if jobs > 0 else default_jobs()


Region = tuple[str, int, int | None]   # contig, start, end (None = to the end)


def has_index(path: str | Path) -> bool:
    p = str(path)
    return os.path.exists(p + ".tbi") or os.path.exists(p + ".csi")


def vcf_regions(path: str | Path, size: int | None = None) -> list[Region]:
    """Regions covering every indexed contig, in FILE order (a tabix index
    lists its contigs in order of appearance, which is what the stages
    restitch — header order is not guaranteed to be file order). The last
    region of a contig is open-ended so records past the declared length
    are kept."""
    size = size or REGION_SIZE
    with pysam.VariantFile(str(path)) as vf:
        lengths = {c.name: c.length or 0 for c in vf.header.contigs.values()}
        names = list(vf.index)
    regions: list[Region] = []
    for name in names:
        length = lengths.get(name, 0)
        start = 0
        while start + size < length:
            regions.append((name, start, start + size))
            start += size
        regions.append((name, start, None))
    return regions


def split_regions(path: str | Path, jobs: int) -> list[Region | None]:
    """The regions to map over with ``jobs`` workers, or [None] (the whole
    file, one pass) when serial or when the file has no index — warned, as
    a missing .tbi silently cost a 2670-sample PED step 1 h on one core."""
    if jobs <= 1:
        return [None]
    if not has_index(path):
        logger.warning("%s has no .tbi/.csi index: one core instead of %d", path, jobs)
        return [None]
    return vcf_regions(path)


def fetch_region(vf: pysam.VariantFile, region: Region) -> Iterator:
    """Records whose start lies in the region, in file order."""
    contig, start, end = region
    try:
        it = vf.fetch(contig, start, end) if end is not None else vf.fetch(contig, start)
    except ValueError:  # contig absent from the index: no record
        return
    for rec in it:
        if rec.start < start or (end is not None and rec.start >= end):
            continue
        yield rec


def pmap(fn: Callable, items: Iterable, jobs: int) -> Iterator:
    """``map`` over a fork-context process pool, results in input order;
    serial in-process when ``jobs <= 1``."""
    items = list(items)
    if jobs <= 1 or len(items) <= 1:
        yield from map(fn, items)
        return
    with ProcessPoolExecutor(
        max_workers=min(jobs, len(items)), mp_context=get_context("fork")
    ) as pool:
        yield from pool.map(fn, items)


# the empty block every BGZF file ends with (SAM spec, 4.1.2)
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def concat_bgzf(parts: list[str | Path], out: str | Path) -> None:
    """Byte-level concatenation of BGZF files into one (a BGZF stream is a
    series of independent blocks: the header part and each worker's records
    are joined without recompression); interior EOF blocks are dropped."""
    with open(out, "wb") as fo:
        for p in parts:
            data = Path(p).read_bytes()
            if data.endswith(BGZF_EOF):
                data = data[: -len(BGZF_EOF)]
            fo.write(data)
        fo.write(BGZF_EOF)


def concat_text(parts: list[str | Path], out: str | Path) -> None:
    with open(out, "wb") as fo:
        for p in parts:
            with open(p, "rb") as fi:
                while chunk := fi.read(1 << 24):
                    fo.write(chunk)
