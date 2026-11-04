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

SEMANTICS: PoN mode reports deviation from the population norm (gCNV-style).
Common CNVs (carrier fraction high enough to move the median, or polymorphic
enough to trip the MAD mask) leave the DEPTH channel — they stay detectable
through PE/SR junctions when a CRAM is scanned. Use --pon for rare-variant
analyses (de novo, case-control); build frequency catalogues without it.
Validated on GHFC chr20: 20-founder panel cut depth segments 264 -> 14, with
12/14 survivors rare (<=1/20 founders) and known common DELs absorbed.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np
import zarr

from ..extract.store import EvidenceStore, _create_array
from ..vcfout import span, spec_coords

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
    denoise_k: int = 0,
    max_svd_samples: int = 2000,
) -> Path:
    """Build the PoN; ``denoise_k > 0`` adds SVD components (SavvyCNV/XHMM style).

    The SVD captures the top-k axes of residual systematic variation of
    log2(sample ratio / panel median) across the panel; at call time the
    sample's projection onto those axes is subtracted. Real CNVs are sparse
    and project weakly onto a small k — keep k modest (5-20). Autosomes only.
    """
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

    if denoise_k:
        rng = np.random.default_rng(17)
        svd_subset = (
            stores if len(stores) <= max_svd_samples
            else [stores[i] for i in rng.choice(len(stores), max_svd_samples, replace=False)]
        )
        k = min(denoise_k, max(len(svd_subset) // 5, 1))
        if k < denoise_k:
            logger.warning("PoN SVD: k reduced to %d (only %d samples)", k, len(svd_subset))
        root.attrs["svd_k"] = k

    for contig in contigs:
        n_bins = int(stores[0].root[contig]["depth_raw"].shape[0])
        g = root.create_group(contig)
        is_sex_contig = contig in ("chrX", "chrY", "X", "Y")
        groups = (
            [(sex, subset) for sex, subset in by_sex.items() if subset]
            if is_sex_contig else [("all", stores)]
        )
        median = None
        for label, subset in groups:
            if len(subset) < (min_samples if label == "all" else max(min_samples // 2, 3)):
                logger.warning("PoN %s/%s: only %d samples, skipping", contig, label, len(subset))
                continue
            median, mad = _profile(subset, contig, n_bins)
            _create_array(g, f"median_{label}", median, _BLOCK)
            _create_array(g, f"mad_{label}", mad, _BLOCK)
            g.attrs[f"n_{label}"] = len(subset)

        if denoise_k and not is_sex_contig and median is not None and k >= 1:
            _build_svd(g, svd_subset, contig, median, k)
        logger.info("PoN %s done (%d bins)", contig, n_bins)
    return Path(out_path)


def _build_svd(
    g: zarr.Group,
    subset: list[EvidenceStore],
    contig: str,
    median: np.ndarray,
    k: int,
    min_median: float = 0.25,
) -> None:
    """Top-k right singular vectors of the panel's log-ratio matrix (one contig)."""
    from sklearn.utils.extmath import randomized_svd

    usable = np.flatnonzero(np.isfinite(median) & (median >= min_median))
    if len(usable) < 10 * k:
        logger.warning("PoN SVD %s: too few usable bins, skipping", contig)
        return
    med_u = median[usable]
    matrix = np.empty((len(subset), len(usable)), dtype=np.float32)
    for i, s in enumerate(subset):
        _raw, corr, _mask = s.depth(contig)
        ratio = corr / float(s.meta["mean_depth"])
        matrix[i] = np.log2(np.maximum(ratio[usable], 1e-3) / med_u)
    matrix -= matrix.mean(axis=0, keepdims=True)
    _u, sing, vt = randomized_svd(matrix, n_components=k, random_state=17)
    _create_array(g, "svd_bins", usable.astype(np.uint32), _BLOCK)
    _create_array(g, "svd_components", vt.astype(np.float32), _BLOCK)
    explained = (sing**2) / max((matrix**2).sum(), 1e-9)
    g.attrs["svd_explained"] = [round(float(v), 5) for v in explained]
    logger.info(
        "PoN SVD %s: k=%d on %d bins, explained %s",
        contig, k, len(usable), g.attrs["svd_explained"],
    )


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
        denoise: bool = True,
        level: float = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """(normalized ratio, extended mask).

        Bins are divided by the population median ratio and brought back to
        the sample's baseline ``level`` (ploidy / 2: 0.5 on a male chrX or
        chrY, whose sex-stratified median is already ~0.5 — without it every
        male chrX came out at ratio 1.0 and was called a whole-chromosome
        DUP); bins where the population is too variable (MAD/median > max_cv)
        or systematically aberrant (median < min_median x level) are masked
        as depth-uninformative. When the panel carries SVD components (and
        ``denoise``), the sample's projection onto the panel's
        residual-variation axes is also removed.
        """
        if level <= 0:  # chrY of a female: no baseline, not called
            return ratio, mask
        prof = self.profile(contig, sex)
        if prof is None or len(prof[0]) != len(ratio):
            logger.warning("PoN: no usable profile for %s, ratio left untouched", contig)
            return ratio, mask
        median, mad = prof
        ok = np.isfinite(median) & (median >= min_median * level)
        out = ratio.copy()
        out[ok] = ratio[ok] / median[ok] * level
        bad = ~ok | (np.isfinite(median) & (mad > max_cv * np.maximum(median, 1e-3)))

        if denoise:
            try:
                g = self.root[contig]
                bins = g["svd_bins"][:]
                components = g["svd_components"][:]  # (k, n_usable)
            except KeyError:
                return out, mask | bad
            x = np.log2(np.maximum(out[bins], 1e-3)).astype(np.float32)
            coeffs = components @ x
            correction = components.T @ coeffs
            out[bins] = out[bins] / np.exp2(correction)
            logger.debug(
                "PoN SVD %s: removed projection (|coeffs| %s)",
                contig, np.round(np.abs(coeffs), 2),
            )
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
        spec = spec_coords(vcf.header)
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
                start, end = span(rec, spec)
                out.write(
                    f"{rec.contig}\t{start}\t{end}\t"
                    f"{rec.info.get('SVTYPE', '.')};AF={af:.4f};{reason}\n"
                )
                n += 1
    logger.info("artifact catalogue: %d loci -> %s", n, out_bed)
    return n
