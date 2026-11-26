"""`svengers filter`: annotate a single-sample VCF with context, hard-filter
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
    ("MODEL_VER", "1", "String", "SVengers filter model version"),
    ("REFMINOR", "1", "Float",
     "Fraction of the call's depth bins where GRCh38 carries the minor allele (most of the panel "
     "of normals deleted there): read against the reference, a deletion most people carry"),
    ("SQ_OVERRIDE", "1", "String",
     "Reason the calibrated SQ was overridden (LARGE_DEL_DEPTH: deletion beyond the "
     "model's training size range whose chromosome-level depth confirms a copy loss; "
     "LARGE_DUP_BAF: large duplication whose heterozygous sites show the allelic "
     "imbalance of a copy gain; TERMINAL_DEL_DEPTH: deletion reaching a contig end, "
     "with a single breakpoint and so little junction evidence, whose depth shows a clean "
     "copy loss; UNIQUE_DEL_DEPTH: depth-only deletion >= 10 kb outside repeats whose depth "
     "shows a clean copy loss, in a genome without coverage waves; JUNCTION_DUP_DEPTH: "
     "duplication of 1-100 kb seen by both depth and junction reads, CN >= 3, in a genome "
     "without strong coverage waves; UNIQUE_DUP_DEPTH: duplication >= 10 kb outside repeats "
     "whose depth shows a clean copy gain, in a genome without strong coverage waves)"),
]


def reaches_contig_end(store, contig: str, start: int, end: int, length: int | None) -> bool:
    """A loss reaching a contig end, as `bnd` reads it: within TERMINAL_BP of
    it, with at most TERMINAL_TAIL_BP of callable (unmasked) sequence past it."""
    import numpy as np

    from ..bnd import TERMINAL_BP, TERMINAL_TAIL_BP

    if not length or not store.has_contig(contig):
        return False
    qter = end >= length - TERMINAL_BP
    if not (qter or start <= TERMINAL_BP):
        return False
    _, _, mask = store.depth(contig)
    pos = np.nonzero(~np.asarray(mask, dtype=bool))[0] * store.bin_size
    tail = (pos >= end).sum() if qter else (pos + store.bin_size <= start).sum()
    return tail * store.bin_size <= TERMINAL_TAIL_BP


def sq_override(svtype: str, size: int, feats: dict, cf, terminal: bool = False) -> str | None:
    """Why a low SQ does not reject this call, or None (see FilterConfig).
    `terminal`: the call reaches a contig end (`reaches_contig_end`)."""
    if (svtype == "DEL" and size >= cf.large_del_rescue_size
            and feats["fc_chrom"] < cf.large_del_rescue_fc_chrom
            and (feats["maskfrac"] < cf.large_del_rescue_max_mask
                 or feats["n_bins_informative"] >= cf.informative_core_bins)):
        return "LARGE_DEL_DEPTH"
    if (svtype == "DEL" and terminal and size >= cf.terminal_del_rescue_size
            and 0 <= feats["cn"] <= 1
            and feats["dhffc"] <= cf.terminal_del_rescue_max_ratio
            and feats["frac_consistent"] >= cf.terminal_del_rescue_min_consistent
            and feats["maskfrac"] < cf.terminal_del_rescue_max_mask):
        return "TERMINAL_DEL_DEPTH"
    if (svtype == "DEL" and size >= cf.unique_del_rescue_size
            and feats["src_both"] == 0 and feats["src_junction"] == 0
            and 0 <= feats["cn"] <= 1
            and feats["dhffc"] <= cf.unique_del_rescue_max_ratio
            and feats["frac_consistent"] >= cf.unique_del_rescue_min_consistent
            and feats["segdup_ovl"] < cf.unique_del_rescue_max_repeat
            and feats["lcr_ovl"] < cf.unique_del_rescue_max_repeat
            and feats["maskfrac"] < cf.unique_del_rescue_max_mask
            and 0 <= feats["wave_score"] < cf.unique_del_rescue_max_wave):
        return "UNIQUE_DEL_DEPTH"
    if (svtype == "DUP" and cf.junction_dup_rescue_min_size <= size < cf.junction_dup_rescue_max_size
            and feats.get("src_both", 0) == 1 and feats["cn"] >= cf.junction_dup_rescue_min_cn
            and 0 <= feats.get("wave_score", -1.0) < cf.junction_dup_rescue_max_wave):
        return "JUNCTION_DUP_DEPTH"
    if (svtype == "DUP" and size >= cf.large_dup_rescue_size and feats["cn"] >= 3
            and feats["baf_n"] >= cf.large_dup_rescue_min_baf_n
            and feats["baf_dev"] >= cf.large_dup_rescue_min_baf_dev
            and feats["segdup_ovl"] < cf.large_dup_rescue_max_segdup
            and feats["maskfrac"] < cf.large_dup_rescue_max_mask):
        return "LARGE_DUP_BAF"
    nan = float("nan")
    if (svtype == "DUP" and size >= cf.unique_dup_rescue_size
            and feats.get("dhffc", nan) >= cf.unique_dup_rescue_min_ratio
            and feats.get("frac_consistent", nan) >= cf.unique_dup_rescue_min_consistent
            and feats.get("segdup_ovl", nan) < cf.unique_dup_rescue_max_repeat
            and feats.get("lcr_ovl", nan) < cf.unique_dup_rescue_max_repeat
            and feats.get("maskfrac", nan) < cf.unique_dup_rescue_max_mask
            and 0 <= feats.get("wave_score", -1.0) < cf.unique_dup_rescue_max_wave):
        return "UNIQUE_DUP_DEPTH"
    return None


def genome_overrides(svtypes, spans, feats_all, sq_all, terminal, cf) -> list[str | None]:
    """`sq_override` of every low-SQ call of one genome. UNIQUE_DEL_DEPTH and
    UNIQUE_DUP_DEPTH are dropped for the whole genome when they would fire
    more than `unique_del_rescue_max_per_genome` / `unique_dup_rescue_max_per_genome`
    times: a burst is genome-wide noise."""
    out = [
        sq_override(t, e - s, f, cf, term) if sq < cf.low_sq_threshold else None
        for t, (s, e), f, sq, term in zip(svtypes, spans, feats_all, sq_all, terminal, strict=True)
    ]
    for rule, cap in (("UNIQUE_DEL_DEPTH", cf.unique_del_rescue_max_per_genome),
                      ("UNIQUE_DUP_DEPTH", cf.unique_dup_rescue_max_per_genome)):
        n = out.count(rule)
        if n > cap:
            logger.info("filter: %d %s candidates > %d, rule off for this genome", n, rule, cap)
            out = [None if o == rule else o for o in out]
    return out


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
    bridged = None
    if config.segment.bridge_max_gap_bins > 0:
        # calls made before call-time bridging (v1.12.0) still carry one event
        # in abutting depth pieces: join them as `call` would (genotype/bridge.py)
        from ..genotype.bridge import bridge_calls

        out_name = Path(vcf_out).name
        bridged = Path(vcf_out).with_name(out_name.split(".vcf")[0] + ".bridged.tmp.vcf.gz")
        if bridge_calls(vcf_in, store_path, bridged, config)["joined"]:
            vcf_in = bridged
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
    terminal: list[bool] = []
    with pysam.VariantFile(str(vcf_in)) as vcf:
        sample = list(vcf.header.samples)[0]
        lengths = {c: vcf.header.contigs[c].length for c in vcf.header.contigs}
        cf = config.filter
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
            terminal.append(
                svtype == "DEL" and end - start >= cf.terminal_del_rescue_size
                and reaches_contig_end(store, rec.contig, start, end, lengths.get(rec.contig))
            )

    rm_cache: dict = {}

    def refminor_frac(contig: str, start: int, end: int) -> float:
        if contig not in rm_cache:
            rm_cache[contig] = store.pon_refminor(contig) if store.has_contig(contig) else None
        a = rm_cache[contig]
        if a is None:
            return 0.0
        b0, b1 = start // store.bin_size, max(-(-end // store.bin_size), start // store.bin_size + 1)
        seg = a[b0:b1]
        return float(seg.mean()) if len(seg) else 0.0

    sq_all = model.score_batch(feats_all, svtypes) if model is not None else None
    overrides = (genome_overrides(svtypes, spans, feats_all, sq_all, terminal, cf)
                 if sq_all is not None else [None] * len(feats_all))

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
                rm = refminor_frac(rec.contig, *spans[i])
                if rm > 0:
                    new_rec.info["REFMINOR"] = round(rm, 3)

                if sq_all is not None:
                    sq = float(sq_all[i])
                    new_rec.samples[sample]["SQ"] = round(sq, 2)
                    new_rec.info["MODEL_VER"] = model.version
                    override = overrides[i]
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
    if bridged is not None:
        for p in (bridged, Path(f"{bridged}.tbi")):
            p.unlink(missing_ok=True)
    logger.info(
        "filter: %(total)d calls, %(pass)d PASS, %(hard_filtered)d hard-filtered, "
        "%(low_sq)d LowSQ, %(sq_override)d large-event SQ overrides", counts,
    )
    return counts
