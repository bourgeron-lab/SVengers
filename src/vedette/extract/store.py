"""The per-sample evidence store.

One Zarr group per sample holding everything downstream stages need, so that
filtering and cohort re-genotyping NEVER re-read the CRAM:

    /                      attrs: schema, sample, reference_id, bin_size,
                                  insert_mu, insert_sigma, mean_depth, ...
    <contig>/depth_raw     float32[n_bins]   raw binned depth
    <contig>/depth_corr    float32[n_bins]   GC-corrected depth (same units)
    <contig>/mask          bool[n_bins]      bins excluded from depth modeling
    <contig>/baf/{pos,vaf,dp}                B-allele fraction at het-callable sites
    <contig>/pe/{pos1,pos2,mate_tid,orient,mapq}   discordant pairs, sorted by pos1
    <contig>/sr/{pos,side,mapq,cliplen}            clipped reads, sorted by pos

PE orientation codes: 0 = DEL-like (+/- long insert), 1 = DUP-like everted (-/+),
2 = INV-like (same strand), 3 = inter-contig. SR side: 0 = clip on the left of
pos (junction ends here), 1 = clip on the right of pos (junction starts here).

Estimated size: 80-120 MB per 35x genome at 100 bp bins -> ~1 TB for 10k samples.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import zarr

SCHEMA_VERSION = "0.1"

PE_DEL, PE_DUP, PE_INV, PE_INTER = 0, 1, 2, 3
SR_LEFT_CLIP, SR_RIGHT_CLIP = 0, 1

_DEPTH_CHUNK = 262_144   # bins per chunk (~26 Mb of genome at 100 bp bins)
_EVENT_CHUNK = 1_048_576


def _create_array(group: zarr.Group, name: str, data: np.ndarray, chunk: int) -> None:
    """Zarr v2/v3 compatible array creation."""
    chunks = (min(chunk, max(len(data), 1)),)
    if hasattr(group, "create_array"):  # zarr >= 3
        arr = group.create_array(name, shape=data.shape, dtype=data.dtype, chunks=chunks)
        if len(data):
            arr[:] = data
    else:  # zarr 2.x
        group.create_dataset(name, data=data, chunks=chunks)


@dataclass
class InsertModel:
    mu: float
    sigma: float
    n_pairs: int
    read_length: float = 150.0

    def discordant_threshold(self, z: float) -> float:
        return self.mu + z * self.sigma


class EvidenceStoreWriter:
    def __init__(self, path: str | Path, sample: str, reference_id: str, bin_size: int):
        self.path = Path(path)
        self.root = zarr.open_group(str(self.path), mode="w")
        self.root.attrs.update(
            schema_version=SCHEMA_VERSION,
            sample=sample,
            reference_id=reference_id,
            bin_size=bin_size,
            created=time.strftime("%Y-%m-%dT%H:%M:%S"),
        )

    def set_meta(self, **kwargs: Any) -> None:
        self.root.attrs.update(**kwargs)

    def _contig(self, contig: str) -> zarr.Group:
        try:
            return self.root[contig]
        except KeyError:
            return self.root.create_group(contig)

    def write_depth(
        self, contig: str, raw: np.ndarray, corrected: np.ndarray, mask: np.ndarray
    ) -> None:
        g = self._contig(contig)
        _create_array(g, "depth_raw", raw.astype(np.float32), _DEPTH_CHUNK)
        _create_array(g, "depth_corr", corrected.astype(np.float32), _DEPTH_CHUNK)
        _create_array(g, "mask", mask.astype(bool), _DEPTH_CHUNK)

    def write_baf(self, contig: str, pos: np.ndarray, vaf: np.ndarray, dp: np.ndarray) -> None:
        order = np.argsort(pos, kind="stable")
        g = self._contig(contig).create_group("baf")
        _create_array(g, "pos", pos[order].astype(np.uint32), _EVENT_CHUNK)
        _create_array(g, "vaf", vaf[order].astype(np.float32), _EVENT_CHUNK)
        _create_array(g, "dp", dp[order].astype(np.uint16), _EVENT_CHUNK)

    def write_pe(
        self,
        contig: str,
        pos1: np.ndarray,
        pos2: np.ndarray,
        mate_tid: np.ndarray,
        orient: np.ndarray,
        mapq: np.ndarray,
    ) -> None:
        order = np.argsort(pos1, kind="stable")
        g = self._contig(contig).create_group("pe")
        _create_array(g, "pos1", pos1[order].astype(np.uint32), _EVENT_CHUNK)
        _create_array(g, "pos2", pos2[order].astype(np.uint32), _EVENT_CHUNK)
        _create_array(g, "mate_tid", mate_tid[order].astype(np.int32), _EVENT_CHUNK)
        _create_array(g, "orient", orient[order].astype(np.uint8), _EVENT_CHUNK)
        _create_array(g, "mapq", mapq[order].astype(np.uint8), _EVENT_CHUNK)

    def write_sr(
        self, contig: str, pos: np.ndarray, side: np.ndarray, mapq: np.ndarray, cliplen: np.ndarray
    ) -> None:
        order = np.argsort(pos, kind="stable")
        g = self._contig(contig).create_group("sr")
        _create_array(g, "pos", pos[order].astype(np.uint32), _EVENT_CHUNK)
        _create_array(g, "side", side[order].astype(np.uint8), _EVENT_CHUNK)
        _create_array(g, "mapq", mapq[order].astype(np.uint8), _EVENT_CHUNK)
        _create_array(g, "cliplen", cliplen[order].astype(np.uint16), _EVENT_CHUNK)


class EvidenceStore:
    """Read-side API — the contract used by filter and cohort stages."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.root = zarr.open_group(str(self.path), mode="r")
        self.meta: dict[str, Any] = dict(self.root.attrs)
        if self.meta.get("schema_version") != SCHEMA_VERSION:
            raise ValueError(
                f"{path}: evidence store schema "
                f"{self.meta.get('schema_version')!r} != expected {SCHEMA_VERSION!r}"
            )
        self.bin_size: int = int(self.meta["bin_size"])

    @property
    def sample(self) -> str:
        return self.meta["sample"]

    @property
    def contigs(self) -> list[str]:
        return sorted(self.root.group_keys()) if hasattr(self.root, "group_keys") else sorted(self.root)

    def insert_model(self) -> InsertModel:
        return InsertModel(
            mu=float(self.meta["insert_mu"]),
            sigma=float(self.meta["insert_sigma"]),
            n_pairs=int(self.meta.get("insert_n_pairs", 0)),
        )

    def has_contig(self, contig: str) -> bool:
        try:
            self.root[contig]
            return True
        except KeyError:
            return False

    def depth(
        self, contig: str, start: int | None = None, end: int | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(raw, corrected, mask) for bins overlapping [start, end) — whole contig if None."""
        g = self.root[contig]
        if start is None:
            sl = slice(None)
        else:
            sl = slice(start // self.bin_size, -(-(end) // self.bin_size))
        return g["depth_raw"][sl], g["depth_corr"][sl], g["mask"][sl]

    def _ranged(self, group: zarr.Group, pos_name: str, start: int, end: int) -> slice:
        pos = group[pos_name][:]
        return slice(
            int(np.searchsorted(pos, start, side="left")),
            int(np.searchsorted(pos, end, side="right")),
        )

    def has_junction_evidence(self, contig: str) -> bool:
        """True if this store recorded PE/SR arrays for the contig (CRAM scan ran)."""
        try:
            self.root[contig]["pe"]
            return True
        except KeyError:
            return False

    def baf(self, contig: str, start: int, end: int) -> dict[str, np.ndarray]:
        g = self.root[contig]
        try:
            baf = g["baf"]
        except KeyError:
            return {"pos": np.array([], np.uint32), "vaf": np.array([], np.float32), "dp": np.array([], np.uint16)}
        sl = self._ranged(baf, "pos", start, end)
        return {name: baf[name][sl] for name in ("pos", "vaf", "dp")}

    def pe(self, contig: str, start: int, end: int) -> dict[str, np.ndarray]:
        """Discordant pairs whose LEFT anchor falls in [start, end)."""
        g = self.root[contig]
        try:
            pe = g["pe"]
        except KeyError:
            return {k: np.array([]) for k in ("pos1", "pos2", "mate_tid", "orient", "mapq")}
        sl = self._ranged(pe, "pos1", start, end)
        return {name: pe[name][sl] for name in ("pos1", "pos2", "mate_tid", "orient", "mapq")}

    def sr(self, contig: str, start: int, end: int) -> dict[str, np.ndarray]:
        g = self.root[contig]
        try:
            sr = g["sr"]
        except KeyError:
            return {k: np.array([]) for k in ("pos", "side", "mapq", "cliplen")}
        sl = self._ranged(sr, "pos", start, end)
        return {name: sr[name][sl] for name in ("pos", "side", "mapq", "cliplen")}
