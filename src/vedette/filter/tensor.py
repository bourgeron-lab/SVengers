"""Per-call tensors for the CNN second opinion (ROADMAP 1sexdecies).

The filter model sees summaries of each call (``dhffc``, ``depth_cv``…), never
the shape of the signal: a clean two-step deletion (two drops to 0.5, flat
flanks) and a GC wave can share the same summaries — the 9 kb SHANK3 deletion
of PMS, obvious in IGV, scored SQ 0.38. This turns each call into a
fixed-length multi-channel track a 1D CNN can read. Everything comes from the
evidence store (the invariant: never the CRAM); the store already holds the
signal condensed, so nothing is rasterized to images.

Window: the call plus one call length of flank on each side (at least
``MIN_FLANK``), resampled to ``L`` points. Channels, in ``CHANNELS`` order:

- ``depth``: log2 copy ratio against the contig's normal ploidy, mean of the
  unmasked 1 kb bins under the point (0 where none);
- ``depth_masked``: fraction of masked or off-contig bins;
- ``pon_median``, ``pon_mad``: the panel of normals the store was called
  with (300 unrelated founders: log2 median ratio and MAD per bin, the
  sex-matched profile on chrX/chrY) — the cohort context. A dip the whole
  cohort shares, or bins where it scatters, are the systematic artifacts a
  single genome cannot tell from a CNV. 0 without a panel;
- ``het_density``, ``baf_dev``: heterozygous sites (VAF 0.1-0.9) per kb and
  their mean |VAF - 0.5| (a one-copy loss has no het sites, a gain shifts
  them to 1/3-2/3);
- ``pe_del``, ``pe_dup``: ends of DEL-like and DUP-like discordant pairs
  (both mates), ``sr_left``, ``sr_right``: soft clips by side, ``cig_del``:
  large CIGAR deletions (both ends) — all per point, scaled to 30x;
- ``ins_z``: local insert-size shift, (mean - mu) / sigma;
- ``in_call``: fraction of the point inside the call.

Scalars (``SCALARS``) go next to the convolution output: size, type, source,
contig ploidy, the context annotations of the filtered VCF, wave score,
coverage. Like the feature schema, ``CHANNELS`` and ``SCALARS`` are a model
contract: append only and bump ``TENSOR_VERSION``.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import pysam

from ..extract.store import PE_DEL, PE_DUP, SR_LEFT_CLIP, SR_RIGHT_CLIP, EvidenceStore

logger = logging.getLogger(__name__)

TENSOR_VERSION = 1
L = 256
MIN_FLANK = 5_000
MIN_MAPQ = 10                # as config.min_mapq_evidence
INS_BIN = 100                # as config.insert_bin_size
DEPTH_LOG2_RANGE = (-3.0, 2.0)
INS_Z_RANGE = (-5.0, 5.0)
SCALE_DEPTH = 30.0           # event counts are scaled to this coverage

CHANNELS: list[str] = [
    "depth", "depth_masked", "pon_median", "pon_mad",
    "het_density", "baf_dev",
    "pe_del", "pe_dup", "sr_left", "sr_right", "cig_del", "ins_z",
    "in_call",
]
SCALARS: list[str] = [
    "size_log10", "is_dup", "src_depth", "src_both", "src_junction", "src_insert",
    "ploidy", "popcov_sd", "segdup_ovl", "lcr_ovl", "wave_score", "coverage",
]


def window(start: int, end: int) -> tuple[int, int]:
    """[w0, w1): the call plus one call length (>= MIN_FLANK) on each side."""
    flank = max(end - start, MIN_FLANK)
    return start - flank, end + flank


def _cum(values: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(values, dtype=np.float64)])


def _bin_ranges(edges: np.ndarray, bin_size: int, n_bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin index range [lo, hi) under each point, at least one bin per point,
    clipped to the contig, and the unclipped bin count (off-contig bins
    count as missing)."""
    lo = np.floor_divide(edges[:-1], bin_size)
    hi = np.maximum(-(-edges[1:] // bin_size), lo + 1)
    return np.clip(lo, 0, n_bins), np.clip(hi, 0, n_bins), hi - lo


class _ContigTracks:
    """One contig of one store, loaded once and cut into windows."""

    def __init__(self, store: EvidenceStore, contig: str, pon_group, normal: int):
        self.bin_size = store.bin_size
        self.normal = normal
        _raw, corr, mask = store.depth(contig)
        corr = np.asarray(corr, dtype=np.float64)
        valid = ~np.asarray(mask, dtype=bool) & (corr > 0)
        mean_depth = float(store.meta["mean_depth"])
        # copy ratio against the contig's normal ploidy (male chrX: 0.5 -> 0)
        ratio = np.where(valid, corr / mean_depth / (normal / 2.0), 1.0)
        log2 = np.clip(np.log2(np.maximum(ratio, 1e-3)), *DEPTH_LOG2_RANGE)
        self.n_bins = len(corr)
        self.c_log2, self.c_valid = _cum(np.where(valid, log2, 0.0)), _cum(valid)

        self.c_pon_med = self.c_pon_mad = None
        if pon_group is not None:
            med, mad = pon_group
            med = np.asarray(med[: self.n_bins], dtype=np.float64)
            mad = np.asarray(mad[: self.n_bins], dtype=np.float64)
            pon_log2 = np.clip(np.log2(np.maximum(med, 1e-3) / (normal / 2.0)), *DEPTH_LOG2_RANGE)
            self.c_pon_med, self.c_pon_mad = _cum(pon_log2), _cum(mad)

        b = store.baf(contig, 0, 1 << 31)
        vaf = np.asarray(b["vaf"], dtype=np.float64)
        het = (vaf > 0.1) & (vaf < 0.9)
        self.het_pos = np.asarray(b["pos"], dtype=np.int64)[het]
        self.het_dev = np.abs(vaf[het] - 0.5)

        scale = SCALE_DEPTH / max(mean_depth, 1.0)
        pe = store.pe(contig, 0, 1 << 31)
        keep = np.asarray(pe["mapq"]) >= MIN_MAPQ
        orient = np.asarray(pe["orient"])
        pos1, pos2 = np.asarray(pe["pos1"], np.int64), np.asarray(pe["pos2"], np.int64)
        # both mates: a deletion's pairs pile up on each side of it
        self.pe_del, self.pe_dup = (np.concatenate([pos1[sel], pos2[sel]])
                                    for sel in (keep & (orient == PE_DEL), keep & (orient == PE_DUP)))
        sr = store.sr(contig, 0, 1 << 31)
        keep = np.asarray(sr["mapq"]) >= MIN_MAPQ
        side = np.asarray(sr["side"])
        pos = np.asarray(sr["pos"], np.int64)
        self.sr_left, self.sr_right = pos[keep & (side == SR_LEFT_CLIP)], pos[keep & (side == SR_RIGHT_CLIP)]
        cig = store.cig(contig, 0, 1 << 31)
        cpos = np.asarray(cig["pos"], np.int64)
        self.cig_ends = np.concatenate([cpos, cpos + np.asarray(cig["svlen"], np.int64)])
        self.scale = scale

        ins = store.ins(contig)
        self.c_ins_sum = self.c_ins_n = None
        if ins is not None:
            self.c_ins_sum, self.c_ins_n = _cum(ins["sum"]), _cum(ins["n"])
            self.n_ins = len(ins["n"])
        im = store.insert_model()
        self.ins_mu, self.ins_sigma = im.mu, max(im.sigma, 1.0)

    def tensor(self, start: int, end: int) -> np.ndarray:
        w0, w1 = window(start, end)
        edges = np.linspace(w0, w1, L + 1).round().astype(np.int64)
        span_kb = np.diff(edges) / 1000.0
        x = np.zeros((len(CHANNELS), L), dtype=np.float32)
        ch = {name: i for i, name in enumerate(CHANNELS)}

        lo, hi, n_all = _bin_ranges(edges, self.bin_size, self.n_bins)
        n_valid = self.c_valid[hi] - self.c_valid[lo]
        with np.errstate(invalid="ignore", divide="ignore"):
            x[ch["depth"]] = np.where(n_valid > 0, (self.c_log2[hi] - self.c_log2[lo]) / n_valid, 0.0)
            x[ch["depth_masked"]] = 1.0 - n_valid / n_all
            if self.c_pon_med is not None:
                n_in = np.maximum(hi - lo, 1)
                x[ch["pon_median"]] = np.where(hi > lo, (self.c_pon_med[hi] - self.c_pon_med[lo]) / n_in, 0.0)
                x[ch["pon_mad"]] = np.where(hi > lo, (self.c_pon_mad[hi] - self.c_pon_mad[lo]) / n_in, 0.0)

        n_het, _ = np.histogram(self.het_pos, bins=edges)
        dev_sum, _ = np.histogram(self.het_pos, bins=edges, weights=self.het_dev)
        x[ch["het_density"]] = n_het / span_kb
        x[ch["baf_dev"]] = np.where(n_het > 0, dev_sum / np.maximum(n_het, 1), 0.0)

        for name, pos in (("pe_del", self.pe_del), ("pe_dup", self.pe_dup), ("sr_left", self.sr_left),
                          ("sr_right", self.sr_right), ("cig_del", self.cig_ends)):
            x[ch[name]] = np.histogram(pos, bins=edges)[0] * self.scale

        if self.c_ins_sum is not None:
            ilo, ihi, _ = _bin_ranges(edges, INS_BIN, self.n_ins)
            n = self.c_ins_n[ihi] - self.c_ins_n[ilo]
            mean = (self.c_ins_sum[ihi] - self.c_ins_sum[ilo]) / np.maximum(n, 1)
            x[ch["ins_z"]] = np.where(n > 0, np.clip((mean - self.ins_mu) / self.ins_sigma, *INS_Z_RANGE), 0.0)

        inside = np.clip(np.minimum(edges[1:], end) - np.maximum(edges[:-1], start), 0, None)
        x[ch["in_call"]] = inside / np.diff(edges)
        return x


def _pon_profile(pon_root, contig: str, normal: int):
    """(median, mad) arrays of the panel for this contig, sex-matched on chrX/chrY."""
    if pon_root is None:
        return None
    try:
        g = pon_root[contig]
    except KeyError:
        return None
    keys = list(g.array_keys())
    if "median_all" in keys:
        return g["median_all"], g["mad_all"]
    sex = "male" if normal == 1 or contig == "chrY" else "female"
    if f"median_{sex}" in keys:
        return g[f"median_{sex}"], g[f"mad_{sex}"]
    return None


def _info(rec: pysam.VariantRecord, key: str, default):
    # pysam raises on a key the header does not declare (VCFs filtered before
    # the context annotations existed)
    return rec.info.get(key, default) if key in rec.header.info else default


def scalars(rec: pysam.VariantRecord, normal: int, wave: float, coverage: float) -> np.ndarray:
    src = _info(rec, "SOURCE", "DEPTH")
    vals = {
        "size_log10": math.log10(max(rec.stop - rec.start, 1)),
        "is_dup": float(rec.info["SVTYPE"] == "DUP"),
        "src_depth": float(src == "DEPTH"),
        "src_both": float(src == "BOTH"),
        "src_junction": float(src == "JUNCTION"),
        "src_insert": float(src == "INSERT"),
        "ploidy": normal / 2.0,
        "popcov_sd": float(_info(rec, "POPCOV_SD", -1.0)),
        "segdup_ovl": float(_info(rec, "SEGDUP_OVL", 0.0)),
        "lcr_ovl": float(_info(rec, "LCR_OVL", 0.0)),
        "wave_score": wave,
        "coverage": coverage / SCALE_DEPTH,
    }
    return np.array([vals[k] for k in SCALARS], dtype=np.float32)


def run_tensors(
    vcf_in: str | Path,
    store_path: str | Path,
    out: str | Path,
    *,
    min_size: int = 1000,
    filters: tuple[str, ...] | None = ("PASS", "LowSQ"),
) -> int:
    """Tensors of the calls of a filtered single-sample VCF -> compressed NPZ.

    Default scope: calls past the hard filters (PASS or LowSQ), >= 1 kb — the
    calls the model scores and where depth has a shape to read (~1-2k per
    genome). ``filters=None`` keeps every call."""
    import zarr

    store = EvidenceStore(store_path)
    pon = store.meta.get("pon")
    pon_root = zarr.open(pon, mode="r") if pon and Path(pon).exists() else None
    if pon_root is not None and int(pon_root.attrs.get("bin_size", store.bin_size)) != store.bin_size:
        logger.warning("panel bin size differs from the store's: cohort context channels left at 0")
        pon_root = None
    ploidy = store.meta.get("ploidy") or {}
    wave = store.meta.get("wave_score")
    wave = float(wave) if wave is not None else -1.0
    coverage = float(store.meta.get("coverage", store.meta.get("mean_depth", 30.0)))

    xs, ss, meta = [], [], {k: [] for k in ("call_id", "contig", "start", "end", "svtype", "filter", "sq")}
    tracks: dict[str, _ContigTracks] = {}
    with pysam.VariantFile(str(vcf_in)) as vcf:
        sample = list(vcf.header.samples)[0]
        for rec in vcf:
            flt = ",".join(rec.filter.keys()) or "PASS"
            if rec.stop - rec.start < min_size or (filters is not None and flt not in filters):
                continue
            if not store.has_contig(rec.contig):
                continue
            normal = int(ploidy.get(rec.contig, 2))
            if normal < 1:
                continue
            if rec.contig not in tracks:
                tracks.clear()  # one contig in memory at a time (records are sorted)
                tracks[rec.contig] = _ContigTracks(store, rec.contig,
                                                   _pon_profile(pon_root, rec.contig, normal), normal)
            xs.append(tracks[rec.contig].tensor(rec.start, rec.stop))
            ss.append(scalars(rec, normal, wave, coverage))
            sq = rec.samples[sample].get("SQ")
            for k, v in (("call_id", rec.id), ("contig", rec.contig), ("start", rec.start),
                         ("end", rec.stop), ("svtype", rec.info["SVTYPE"]), ("filter", flt),
                         ("sq", np.nan if sq is None else sq)):
                meta[k].append(v)

    n = len(xs)
    np.savez_compressed(
        out,
        x=np.stack(xs).astype(np.float16) if n else np.zeros((0, len(CHANNELS), L), np.float16),
        s=np.stack(ss) if n else np.zeros((0, len(SCALARS)), np.float32),
        channels=np.array(CHANNELS), scalars=np.array(SCALARS),
        tensor_version=TENSOR_VERSION, sample=store.sample, has_pon=pon_root is not None,
        **{k: np.array(v) for k, v in meta.items()},
    )
    return n
