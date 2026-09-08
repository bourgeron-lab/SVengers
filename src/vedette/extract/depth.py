"""Binned read depth: fast path (precomputed bedgraph) and GC correction.

The self-contained CRAM depth counter lives in ``scan.py``; this module handles
the depth *track* itself: loading it from a tabix-indexed bedgraph (dense bins,
e.g. the GHFC ``{ID}.by1000.bedgraph.gz`` files) and correcting it for GC bias.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pysam

logger = logging.getLogger(__name__)


def read_depth_bedgraph(
    path: str, contig: str, n_bins: int, bin_size: int
) -> np.ndarray:
    """Load a dense fixed-bin bedgraph into a float32 array of n_bins.

    Rows must be aligned on ``bin_size`` boundaries (multi-bin rows are
    broadcast). Missing bins stay 0.
    """
    depth = np.zeros(n_bins, dtype=np.float32)
    with pysam.TabixFile(path) as tbx:
        if contig not in tbx.contigs:
            logger.warning("%s: contig %s absent from bedgraph", path, contig)
            return depth
        for row in tbx.fetch(contig, parser=pysam.asTuple()):
            start, end, value = int(row[1]), int(row[2]), float(row[3])
            b0 = start // bin_size
            b1 = min(-(-end // bin_size), n_bins)
            if b0 < n_bins:
                depth[b0:b1] = value
    return depth


@dataclass
class GCModel:
    """Median depth per GC stratum, fit on unmasked autosomal bins."""

    gc_edges: np.ndarray        # stratum edges, e.g. 0, 0.02, ..., 1
    median_by_gc: np.ndarray    # median raw depth per stratum (NaN if unpopulated)
    global_median: float        # median raw depth over all included bins

    def correction(self, gc: np.ndarray) -> np.ndarray:
        """Multiplicative per-bin correction factor (1.0 where the stratum is empty)."""
        idx = np.clip(np.digitize(gc, self.gc_edges) - 1, 0, len(self.median_by_gc) - 1)
        med = self.median_by_gc[idx]
        factor = np.where(
            np.isfinite(med) & (med > 0), self.global_median / np.where(med > 0, med, 1.0), 1.0
        )
        return factor.astype(np.float32)


def fit_gc_model(
    depth_by_contig: dict[str, np.ndarray],
    gc_by_contig: dict[str, np.ndarray],
    mask_by_contig: dict[str, np.ndarray],
    autosomes: list[str],
    gc_step: float = 0.02,
    min_bins_per_stratum: int = 100,
) -> GCModel:
    depths, gcs = [], []
    for contig in autosomes:
        keep = ~mask_by_contig[contig] & (depth_by_contig[contig] > 0)
        depths.append(depth_by_contig[contig][keep])
        gcs.append(gc_by_contig[contig][keep])
    depth = np.concatenate(depths)
    gc = np.concatenate(gcs)
    if len(depth) == 0:
        raise ValueError("No usable bins to fit the GC model")

    edges = np.arange(0.0, 1.0 + gc_step, gc_step)
    idx = np.clip(np.digitize(gc, edges) - 1, 0, len(edges) - 2)
    medians = np.full(len(edges) - 1, np.nan)
    for stratum in range(len(edges) - 1):
        vals = depth[idx == stratum]
        if len(vals) >= min_bins_per_stratum:
            medians[stratum] = np.median(vals)
    model = GCModel(edges, medians, float(np.median(depth)))
    logger.info(
        "GC model: global median depth %.2f, %d/%d strata populated",
        model.global_median, int(np.isfinite(medians).sum()), len(medians),
    )
    return model


def correct_depth(
    depth: np.ndarray, gc: np.ndarray, model: GCModel
) -> np.ndarray:
    return (depth * model.correction(gc)).astype(np.float32)


def sample_depth_stats(
    corr_by_contig: dict[str, np.ndarray],
    mask_by_contig: dict[str, np.ndarray],
    autosomes: list[str],
) -> dict[str, float]:
    """Diploid mean/MAD depth from corrected autosomal bins (model baseline)."""
    vals = np.concatenate(
        [corr_by_contig[c][~mask_by_contig[c]] for c in autosomes]
    )
    vals = vals[vals > 0]
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med)) * 1.4826)
    return {"mean_depth": med, "depth_mad": mad}
