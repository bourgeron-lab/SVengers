"""`vedette filter`: annotate a single-sample VCF with context, hard-filter
tags, features (optional TSV export) and — when a trained model is provided —
the calibrated SQ score.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pysam

from ..config import Config
from ..extract.store import EvidenceStore
from ..vcfout import span, spec_coords
from .annotate import ContextAnnotator
from .features import FeatureExtractor
from .hard_filters import hard_filter_tags
from .model import FilterModel
from .schema import FEATURES

logger = logging.getLogger(__name__)

_NEW_INFO = [
    ("SEGDUP_OVL", "1", "Float", "Fraction of call overlapping segmental duplications"),
    ("LCR_OVL", "1", "Float", "Fraction of call overlapping low-complexity/self-chain regions"),
    ("EXCLUDE_OVL", "1", "Float", "Fraction of call overlapping the exclusion mask"),
    ("POPCOV_SD", "1", "Float", "Mean population coverage stdev over the call"),
    ("MODEL_VER", "1", "String", "Vedette filter model version"),
    ("SQ_OVERRIDE", "1", "String",
     "Reason the calibrated SQ was overridden (LARGE_DEL_DEPTH: deletion beyond the "
     "model's training size range whose chromosome-level depth confirms a copy loss; "
     "LARGE_DUP_BAF: large duplication whose heterozygous sites show the allelic "
     "imbalance of a copy gain)"),
]


def sq_override(svtype: str, size: int, feats: dict, cf) -> str | None:
    """Why a low SQ does not reject this call, or None (see FilterConfig)."""
    if (svtype == "DEL" and size >= cf.large_del_rescue_size
            and feats["fc_chrom"] < cf.large_del_rescue_fc_chrom
            and feats["maskfrac"] < cf.large_del_rescue_max_mask):
        return "LARGE_DEL_DEPTH"
    if (svtype == "DUP" and size >= cf.large_dup_rescue_size and feats["cn"] >= 3
            and feats["baf_n"] >= cf.large_dup_rescue_min_baf_n
            and feats["baf_dev"] >= cf.large_dup_rescue_min_baf_dev
            and feats["segdup_ovl"] < cf.large_dup_rescue_max_segdup
            and feats["maskfrac"] < cf.large_dup_rescue_max_mask):
        return "LARGE_DUP_BAF"
    return None


def run_filter(
    vcf_in: str | Path,
    store_path: str | Path,
    vcf_out: str | Path,
    config: Config,
    model_dir: str | Path | None = None,
    features_out: str | Path | None = None,
    artifact_bed: str | Path | None = None,
) -> dict[str, int]:
    from .annotate import IntervalSet

    artifacts = IntervalSet.from_bed(artifact_bed) if artifact_bed else None
    store = EvidenceStore(store_path)
    annotator = ContextAnnotator.from_paths(
        segdup_bed=config.reference.segdup_bed,
        lcr_bed=config.reference.lcr_bed,
        exclude_bed=config.reference.exclude_bed,
        pop_cov_stdev_bw=config.reference.pop_cov_stdev_bw,
    )
    extractor = FeatureExtractor(store, annotator)
    model = FilterModel.load(model_dir) if model_dir else None

    counts = {"total": 0, "pass": 0, "hard_filtered": 0, "low_sq": 0, "sq_override": 0}
    feat_writer = None
    feat_fh = None
    if features_out:
        feat_fh = open(features_out, "w", newline="")
        feat_writer = csv.writer(feat_fh, delimiter="\t")
        feat_writer.writerow(["call_id", "contig", "start", "end", "svtype", *FEATURES])

    # pass 1: extract features and hard-filter tags for every call (the model,
    # when present, scores everything in one vectorized pass afterwards —
    # per-call predict_proba was the runtime bottleneck with deep forests)
    feats_all: list[dict] = []
    tags_all: list[tuple[list[str], list[str]]] = []  # (hard_tags, all_tags)
    svtypes: list[str] = []
    spans: list[tuple[int, int]] = []
    with pysam.VariantFile(str(vcf_in)) as vcf:
        sample = list(vcf.header.samples)[0]
        spec = spec_coords(vcf.header)
        for rec in vcf:
            svtype = rec.info["SVTYPE"]
            start, end = span(rec, spec)
            spans.append((start, end))
            fmt = rec.samples[sample]
            feats = extractor.extract(
                rec.contig,
                start,
                end,
                svtype,
                cn=int(fmt.get("CN") or -1),
                cnq=float(fmt.get("CNQ") or float("nan")),
                gq=int(fmt.get("GQ") or 0),
                pe=int(fmt.get("PE") or 0),
                sr=int(fmt.get("SR") or 0),
                cig=int(fmt.get("CG") or 0),
                ins_z=float(fmt.get("IZ") or 0.0),
                ins_pairs=int(fmt.get("IP") or 0),
                source=rec.info.get("SOURCE", "DEPTH"),
                precise="IMPRECISE" not in rec.info,
                homlen=int(rec.info.get("HOMLEN", -1)),
                baf_n=int(fmt.get("BAFN") or 0),
                baf_dev=float(fmt.get("BAFDEV") or float("nan")),
            )
            tags = hard_filter_tags(feats, end - start, config.filter)
            if (
                artifacts is not None
                and artifacts.overlap_fraction(rec.contig, start, end) > 0.5
            ):
                tags.append("CohortArtifact")
            existing = [f for f in rec.filter.keys() if f != "PASS"]
            feats_all.append(feats)
            tags_all.append((tags, sorted(set(existing) | set(tags))))
            svtypes.append(svtype)

    sq_all = model.score_batch(feats_all, svtypes) if model is not None else None

    # pass 2: re-read the input and write the annotated records in order
    with pysam.VariantFile(str(vcf_in)) as vcf:
        header = vcf.header.copy()
        for name, number, vtype, desc in _NEW_INFO:
            if name not in header.info:
                header.info.add(name, number, vtype, desc)
        sample = list(header.samples)[0]

        out_path = str(vcf_out)
        mode = "wz" if out_path.endswith(".gz") else "w"
        with pysam.VariantFile(out_path, mode, header=header) as out:
            for i, rec in enumerate(vcf):
                counts["total"] += 1
                feats = feats_all[i]
                tags, all_tags = tags_all[i]
                all_tags = list(all_tags)

                new_rec = rec.copy()
                new_rec.translate(header)
                new_rec.info["SEGDUP_OVL"] = round(feats["segdup_ovl"], 4)
                new_rec.info["LCR_OVL"] = round(feats["lcr_ovl"], 4)
                new_rec.info["EXCLUDE_OVL"] = round(feats["exclude_ovl"], 4)
                new_rec.info["POPCOV_SD"] = round(feats["popcov_stdev"], 3)

                if sq_all is not None:
                    sq = float(sq_all[i])
                    new_rec.samples[sample]["SQ"] = round(sq, 2)
                    new_rec.info["MODEL_VER"] = model.version
                    cf = config.filter
                    override = sq_override(svtypes[i], spans[i][1] - spans[i][0], feats, cf)
                    if sq < cf.low_sq_threshold and not override:
                        all_tags.append("LowSQ")
                        counts["low_sq"] += 1
                    elif override and sq < cf.low_sq_threshold:
                        new_rec.info["SQ_OVERRIDE"] = override
                        counts["sq_override"] += 1

                new_rec.filter.clear()
                for tag in all_tags or ["PASS"]:
                    new_rec.filter.add(tag)
                if not all_tags:
                    counts["pass"] += 1
                elif tags:
                    counts["hard_filtered"] += 1
                out.write(new_rec)

                if feat_writer is not None:
                    feat_writer.writerow(
                        [rec.id, rec.contig, *spans[i], svtypes[i]]
                        + [f"{feats[name]:.6g}" for name in FEATURES]
                    )
    if feat_fh is not None:
        feat_fh.close()
    if mode == "wz":
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info(
        "filter: %(total)d calls, %(pass)d PASS, %(hard_filtered)d hard-filtered, "
        "%(low_sq)d LowSQ, %(sq_override)d large-event SQ overrides", counts,
    )
    return counts
