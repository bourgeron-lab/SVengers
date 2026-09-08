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
]


def run_filter(
    vcf_in: str | Path,
    store_path: str | Path,
    vcf_out: str | Path,
    config: Config,
    model_dir: str | Path | None = None,
    features_out: str | Path | None = None,
) -> dict[str, int]:
    store = EvidenceStore(store_path)
    annotator = ContextAnnotator.from_paths(
        segdup_bed=config.reference.segdup_bed,
        lcr_bed=config.reference.lcr_bed,
        exclude_bed=config.reference.exclude_bed,
        pop_cov_stdev_bw=config.reference.pop_cov_stdev_bw,
    )
    extractor = FeatureExtractor(store, annotator)
    model = FilterModel.load(model_dir) if model_dir else None

    counts = {"total": 0, "pass": 0, "hard_filtered": 0, "low_sq": 0}
    feat_writer = None
    feat_fh = None
    if features_out:
        feat_fh = open(features_out, "w", newline="")
        feat_writer = csv.writer(feat_fh, delimiter="\t")
        feat_writer.writerow(["call_id", "contig", "start", "end", "svtype", *FEATURES])

    with pysam.VariantFile(str(vcf_in)) as vcf:
        header = vcf.header.copy()
        for name, number, vtype, desc in _NEW_INFO:
            if name not in header.info:
                header.info.add(name, number, vtype, desc)
        sample = list(header.samples)[0]

        out_path = str(vcf_out)
        mode = "wz" if out_path.endswith(".gz") else "w"
        with pysam.VariantFile(out_path, mode, header=header) as out:
            for rec in vcf:
                counts["total"] += 1
                svtype = rec.info["SVTYPE"]
                fmt = rec.samples[sample]
                feats = extractor.extract(
                    rec.contig,
                    rec.start,
                    rec.stop,
                    svtype,
                    cn=int(fmt.get("CN") or -1),
                    cnq=float(fmt.get("CNQ") or float("nan")),
                    gq=int(fmt.get("GQ") or 0),
                    pe=int(fmt.get("PE") or 0),
                    sr=int(fmt.get("SR") or 0),
                    source=rec.info.get("SOURCE", "DEPTH"),
                    precise="IMPRECISE" not in rec.info,
                    homlen=int(rec.info.get("HOMLEN", -1)),
                    baf_n=int(fmt.get("BAFN") or 0),
                    baf_dev=float(fmt.get("BAFDEV") or float("nan")),
                )

                new_rec = rec.copy()
                new_rec.translate(header)
                new_rec.info["SEGDUP_OVL"] = round(feats["segdup_ovl"], 4)
                new_rec.info["LCR_OVL"] = round(feats["lcr_ovl"], 4)
                new_rec.info["EXCLUDE_OVL"] = round(feats["exclude_ovl"], 4)
                new_rec.info["POPCOV_SD"] = round(feats["popcov_stdev"], 3)

                tags = hard_filter_tags(feats, rec.stop - rec.start, config.filter)
                existing = [f for f in rec.filter.keys() if f != "PASS"]
                all_tags = sorted(set(existing) | set(tags))

                if model is not None:
                    sq = model.score(feats, svtype)
                    new_rec.samples[sample]["SQ"] = round(sq, 2)
                    new_rec.info["MODEL_VER"] = model.version
                    if sq < config.filter.low_sq_threshold:
                        all_tags.append("LowSQ")
                        counts["low_sq"] += 1

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
                        [rec.id, rec.contig, rec.start, rec.stop, svtype]
                        + [f"{feats[name]:.6g}" for name in FEATURES]
                    )
    if feat_fh is not None:
        feat_fh.close()
    if mode == "wz":
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info(
        "filter: %(total)d calls, %(pass)d PASS, %(hard_filtered)d hard-filtered, "
        "%(low_sq)d LowSQ", counts,
    )
    return counts
