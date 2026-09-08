"""Panel of normals: per-bin depth profile and artifact-locus catalogue.

`build_pon` streams the corrected copy-ratio tracks of many evidence stores
(block by block, aligned on the Zarr chunking, so memory stays flat at any
cohort size) into a per-bin median + MAD profile — sex-stratified on chrX/Y.

At call time the PoN divides out residual systematic bias the GC model cannot
see (mappability, reference artifacts, batch-wide oddities) and masks bins
whose population variance makes depth uninformative — the cohort-matched
analogue of ClinSV's popCovStdev track.

`artifacts_bed` derives the recurrent-artifact BED from a cohort VCF
(high founder AF + HWE failure, or nominator/carrier discordance).

IMPORTANT: the panel must be built from UNRELATED individuals (founders), and
enough of them — the median's breakdown point is 50%, so with few or related
samples the panel absorbs real inherited CNVs and normalizes them away
(verified empirically: a trio-only panel erased the trio's own het DELs).
20+ founders is the practical floor; the target sample's relatives must be
excluded.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np
import zarr

from ..extract.store import EvidenceStore, _create_array

logger = logging.getLogger(__name__)

PON_SCHEMA = "0.1"
_BLOCK = 262_144  # bins per streaming block (= depth chunk size)


# ------------------------------------------------------------------------------- build


def _ratio_block(store: EvidenceStore, contig: str, b0: int, b1: int) -> np.ndarray:
    _raw, corr, _mask = store.depth(
        contig, b0 * store.bin_size, b1 * store.bin_size
    )
    return corr / float(store.meta["mean_depth"])


def _sex_of_store(store: EvidenceStore) -> str:
    ploidy = dict(store.meta.get("ploidy", {}))
    return "male" if int(ploidy.get("chrX", 2)) == 1 else "female"


def build_pon(
    store_paths: list[str | Path],
    out_path: str | Path,
    min_samples: int = 8,
) -> Path:
    stores = [EvidenceStore(p) for p in store_paths]
    if not stores:
        raise ValueError("No evidence stores given")
    bin_size = stores[0].bin_size
    reference_id = stores[0].meta.get("reference_id", "")
    for s in stores:
        if s.bin_size != bin_size or s.meta.get("reference_id", "") != reference_id:
            raise ValueError(f"{s.path}: bin size / reference mismatch with the first store")

    root = zarr.open_group(str(out_path), mode="w")
    root.attrs.update(
        schema_version=PON_SCHEMA,
        bin_size=bin_size,
        reference_id=reference_id,
        n_samples=len(stores),
        samples=[s.sample for s in stores],
        created=time.strftime("%Y-%m-%dT%H:%M:%S"),
    )

    contigs = sorted(set.intersection(*(set(s.contigs) for s in stores)))
    by_sex = {"male": [s for s in stores if _sex_of_store(s) == "male"]}
    by_sex["female"] = [s for s in stores if s not in by_sex["male"]]
    logger.info(
        "PoN: %d samples (%d male / %d female), %d contigs, bin %d bp",
        len(stores), len(by_sex["male"]), len(by_sex["female"]), len(contigs), bin_size,
    )
    if len(stores) < 20:
        logger.warning(
            "PoN built from only %d samples: common/inherited CNVs may leak into "
            "the median and erase true calls — use 20+ UNRELATED founders",
            len(stores),
        )

    def _profile(subset: list[EvidenceStore], contig: str, n_bins: int):
        median = np.full(n_bins, np.nan, dtype=np.float32)
        mad = np.full(n_bins, np.nan, dtype=np.float32)
        for b0 in range(0, n_bins, _BLOCK):
            b1 = min(b0 + _BLOCK, n_bins)
            stack = np.empty((len(subset), b1 - b0), dtype=np.float32)
            for k, s in enumerate(subset):
                stack[k] = _ratio_block(s, contig, b0, b1)
            med = np.median(stack, axis=0)
            median[b0:b1] = med
            mad[b0:b1] = np.median(np.abs(stack - med[None, :]), axis=0) * 1.4826
        return median, mad

    for contig in contigs:
        n_bins = int(stores[0].root[contig]["depth_raw"].shape[0])
        g = root.create_group(contig)
        is_sex_contig = contig in ("chrX", "chrY", "X", "Y")
        groups = (
            [(sex, subset) for sex, subset in by_sex.items() if subset]
            if is_sex_contig else [("all", stores)]
        )
        for label, subset in groups:
            if len(subset) < (min_samples if label == "all" else max(min_samples // 2, 3)):
                logger.warning("PoN %s/%s: only %d samples, skipping", contig, label, len(subset))
                continue
            median, mad = _profile(subset, contig, n_bins)
            _create_array(g, f"median_{label}", median, _BLOCK)
            _create_array(g, f"mad_{label}", mad, _BLOCK)
            g.attrs[f"n_{label}"] = len(subset)
        logger.info("PoN %s done (%d bins)", contig, n_bins)
    return Path(out_path)


# ------------------------------------------------------------------------------- apply


class PanelOfNormals:
    def __init__(self, path: str | Path):
        self.root = zarr.open_group(str(path), mode="r")
        self.meta = dict(self.root.attrs)
        if self.meta.get("schema_version") != PON_SCHEMA:
            raise ValueError(f"{path}: unsupported PoN schema")
        self.bin_size = int(self.meta["bin_size"])

    def profile(self, contig: str, sex: str = "all") -> tuple[np.ndarray, np.ndarray] | None:
        try:
            g = self.root[contig]
        except KeyError:
            return None
        for label in (sex, "all", "female", "male"):
            try:
                return g[f"median_{label}"][:], g[f"mad_{label}"][:]
            except KeyError:
                continue
        return None

    def normalize(
        self,
        contig: str,
        ratio: np.ndarray,
        mask: np.ndarray,
        sex: str = "all",
        min_median: float = 0.25,
        max_cv: float = 0.30,
    ) -> tuple[np.ndarray, np.ndarray]:
        """(normalized ratio, extended mask).

        Bins are divided by the population median ratio; bins where the
        population is too variable (MAD/median > max_cv) or systematically
        aberrant (median < min_median) are masked as depth-uninformative.
        """
        prof = self.profile(contig, sex)
        if prof is None or len(prof[0]) != len(ratio):
            logger.warning("PoN: no usable profile for %s, ratio left untouched", contig)
            return ratio, mask
        median, mad = prof
        ok = np.isfinite(median) & (median >= min_median)
        out = ratio.copy()
        out[ok] = ratio[ok] / median[ok]
        bad = ~ok | (np.isfinite(median) & (mad > max_cv * np.maximum(median, 1e-3)))
        return out, mask | bad


# --------------------------------------------------------------------------- artifacts


def artifacts_bed(
    cohort_vcf: str | Path,
    out_bed: str | Path,
    min_af: float = 0.03,
    max_hwe_p: float = 1e-6,
    min_gtr_disc: float = 0.5,
    min_nominators: int = 3,
) -> int:
    """Recurrent-artifact BED from a cohort VCF (vedette cohort-filter output)."""
    import pysam

    n = 0
    with pysam.VariantFile(str(cohort_vcf)) as vcf, open(out_bed, "w") as out:
        for rec in vcf:
            af = rec.info.get("AF", (0.0,))
            af = af[0] if isinstance(af, tuple) else af
            hwe = float(rec.info.get("HWE_P", 1.0))
            disc = float(rec.info.get("GTR_DISC", 0.0))
            nominators = int(rec.info.get("NSAMPLES", 0))
            is_artifact = (af >= min_af and hwe <= max_hwe_p) or (
                nominators >= min_nominators and disc >= min_gtr_disc
            )
            if is_artifact:
                reason = "HWE" if hwe <= max_hwe_p else "GTR_DISC"
                out.write(
                    f"{rec.contig}\t{rec.start}\t{rec.stop}\t"
                    f"{rec.info.get('SVTYPE', '.')};AF={af:.4f};{reason}\n"
                )
                n += 1
    logger.info("artifact catalogue: %d loci -> %s", n, out_bed)
    return n
