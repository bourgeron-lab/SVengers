"""Per-call feature extraction from the VCF record + evidence store."""

from __future__ import annotations

import logging
import math

import numpy as np

from ..extract.store import EvidenceStore
from .annotate import ContextAnnotator
from .schema import FEATURES

logger = logging.getLogger(__name__)

_FLANK_BINS_MIN, _FLANK_BINS_MAX = 3, 100


class FeatureExtractor:
    def __init__(self, store: EvidenceStore, annotator: ContextAnnotator):
        self.store = store
        self.annotator = annotator
        self.coverage = float(store.meta.get("coverage", store.meta.get("mean_depth", 30.0)))
        self.mean_depth = float(store.meta["mean_depth"])
        self._chrom_median: dict[str, float] = {}

    def _chrom_median_depth(self, contig: str) -> float:
        if contig not in self._chrom_median:
            _raw, corr, mask = self.store.depth(contig)
            vals = corr[~mask]
            vals = vals[vals > 0]
            self._chrom_median[contig] = float(np.median(vals)) if len(vals) else self.mean_depth
        return self._chrom_median[contig]

    def extract(
        self,
        contig: str,
        start: int,
        end: int,
        svtype: str,
        *,
        cn: int,
        cnq: float,
        gq: int,
        pe: int,
        sr: int,
        cig: int,
        ins_z: float,
        ins_pairs: int,
        source: str,
        precise: bool,
        baf_n: int,
        baf_dev: float,
        homlen: int = -1,
        ploidy: int = 2,
    ) -> dict[str, float]:
        bs = self.store.bin_size
        b0, b1 = start // bs, max(-(-end // bs), start // bs + 1)
        _raw, corr, mask = self.store.depth(contig, start, end)
        inside = corr[~mask]

        n_bins = b1 - b0
        flank_bins = int(np.clip(n_bins, _FLANK_BINS_MIN, _FLANK_BINS_MAX))
        lf0 = max(b0 - flank_bins, 0)
        _r, lcorr, lmask = self.store.depth(contig, lf0 * bs, b0 * bs) if b0 > lf0 else (None, np.array([]), np.array([], bool))
        _r, rcorr, rmask = self.store.depth(contig, b1 * bs, (b1 + flank_bins) * bs)
        left = lcorr[~lmask] if len(lcorr) else np.array([])
        right = rcorr[~rmask] if len(rcorr) else np.array([])
        flanks = np.concatenate([left, right]) if len(left) or len(right) else np.array([])

        call_med = float(np.median(inside)) if len(inside) else float("nan")
        flank_med = float(np.median(flanks)) if len(flanks) else self._chrom_median_depth(contig)
        chrom_med = self._chrom_median_depth(contig)

        expected_ratio = (cn / 2.0) if cn >= 0 else (0.5 if svtype == "DEL" else 1.5)
        ratio_inside = inside / self.mean_depth if len(inside) else np.array([])

        feats: dict[str, float] = {
            "size_log10": math.log10(max(end - start, 1)),
            "is_dup": 1.0 if svtype == "DUP" else 0.0,
            "precise": 1.0 if precise else 0.0,
            "src_both": 1.0 if source == "BOTH" else 0.0,
            "src_junction": 1.0 if source == "JUNCTION" else 0.0,
            "maskfrac": float(mask.mean()) if len(mask) else 1.0,
            "dhffc": call_med / flank_med if flank_med > 0 else float("nan"),
            "fc_chrom": call_med / chrom_med if chrom_med > 0 else float("nan"),
            "fc_residual": (
                abs(float(ratio_inside.mean()) - expected_ratio)
                if len(ratio_inside) else float("nan")
            ),
            "depth_cv": (
                float(inside.std() / inside.mean())
                if len(inside) > 1 and inside.mean() > 0 else float("nan")
            ),
            "frac_consistent": (
                float((np.abs(ratio_inside - expected_ratio) < 0.25).mean())
                if len(ratio_inside) else 0.0
            ),
            "n_bins_informative": float(len(inside)),
            "flank_asym": (
                abs(float(np.median(left)) - float(np.median(right)))
                / max(flank_med, 1e-3)
                if len(left) and len(right) else float("nan")
            ),
            "pe": float(pe),
            "sr": float(sr),
            "cig": float(cig),
            "ins_z": float(ins_z),
            "ins_pairs": float(ins_pairs),
            "evidence_ratio": (pe + sr + cig) / max(self.coverage, 1.0),
            "homlen": float(homlen),
            "baf_n": float(baf_n),
            "baf_dev": baf_dev if baf_n > 0 and not math.isnan(baf_dev) else -1.0,
            "gq": float(gq),
            "cn": float(cn),
            "cnq": float(cnq) if not math.isnan(cnq) else -1.0,
        }
        feats.update(self.annotator.annotate(contig, start, end))
        missing = set(FEATURES) - set(feats)
        if missing:  # keep the schema contract airtight
            raise RuntimeError(f"Feature extractor missing schema features: {sorted(missing)}")
        return feats
