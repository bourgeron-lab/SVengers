"""Single-sample calling pipeline: evidence -> evidence store + VCF.

Two depth sources share one flow:
- fast path: precomputed dense depth bedgraph (e.g. 1 kb) + per-site VAF track;
- self-contained path: single-pass CRAM scan (100 bp read-start counts).

PE/SR junction evidence always comes from the CRAM scan when a CRAM is given;
without one (fast path only) calling is depth+BAF-only. Candidates are the
reconciliation of CN-HMM segments and PE/SR junctions, genotyped by the
unified RD+PE/SR+BAF likelihood.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from .breakpoints.cluster import make_junctions
from .config import Config
from .extract import baf as baf_mod
from .extract import depth as depth_mod
from .extract.store import EvidenceStoreWriter, InsertModel
from .genotype.likelihood import ContigData, Genotyper
from .karyotype import autosomal_aneuploidies, sex_ploidy
from .merge_candidates import reconcile
from .reference import (
    Contig,
    GenomeBins,
    ReferenceTracks,
    build_bin_mask,
    build_tracks,
    read_bed,
    read_fai,
    reference_id,
)
from .segment import hmm
from .variants import CNVCall
from .vcfout import write_vcf

logger = logging.getLogger(__name__)


def estimate_sex_ploidy(
    ratio_by_contig: dict[str, np.ndarray],
    mask_by_contig: dict[str, np.ndarray],
    contigs: dict[str, Contig],
) -> dict[str, int]:
    """Baseline copy number per contig from median X/Y copy ratio."""
    ploidy = {name: 2 for name in ratio_by_contig}
    x_ratio = y_ratio = None
    for name in ratio_by_contig:
        kind = contigs[name].kind
        if kind in ("X", "Y"):
            vals = ratio_by_contig[name][~mask_by_contig[name]]
            vals = vals[vals > 0]
            med = float(np.median(vals)) if len(vals) else 0.0
            if kind == "X":
                x_ratio = med
            else:
                y_ratio = med
    is_male = (x_ratio is not None and x_ratio < 0.75) or (
        y_ratio is not None and y_ratio > 0.25
    )
    for name in ratio_by_contig:
        kind = contigs[name].kind
        if kind == "X":
            ploidy[name] = 1 if is_male else 2
        elif kind == "Y":
            ploidy[name] = 1 if is_male else 0
    logger.info(
        "Sex inference: X ratio %.2f, Y ratio %.2f -> %s",
        -1 if x_ratio is None else x_ratio,
        -1 if y_ratio is None else y_ratio,
        "male" if is_male else "female",
    )
    return ploidy


def pon_profile(kind: str, ploidy: int) -> tuple[str, float]:
    """(panel profile, level) for one contig. The level is the copy number
    of the PROFILE / 2, not the sample's: ``normalize`` returns the ratio in
    diploid units (ratio / panel median x level), which the HMM compares to
    the contig's own ploidy — an XXX chrX on the female profile reads 1.5
    against ploidy 3, an XYY chrY on the male profile 1.0 against ploidy 2."""
    if kind == "X":
        return ("male", 0.5) if ploidy == 1 else ("female", 1.0)
    if kind == "Y":
        return ("male", 0.5) if ploidy > 0 else ("female", 0.0)
    return "all", 1.0


def sex_karyotype(sample, bin_size, raw, mask, baf, pon_path) -> dict:
    """``classify`` on the in-memory depth/BAF of a call (karyotype.py)."""
    import zarr

    from .karyotype import ArraySource, classify, sex_chromosome_stats

    st = sex_chromosome_stats(
        ArraySource(sample, bin_size, raw, mask, baf), zarr.open(str(pon_path), mode="r")
    )
    return classify(st) if st else {"karyotype": "?", "note": "chrX non lisible"}


def run_call(
    sample: str,
    out_dir: str | Path,
    config: Config,
    resources_dir: str | Path,
    cram: str | None = None,
    depth_bedgraph: str | None = None,
    vaf_bedgraph: str | None = None,
    snv_vcf: str | None = None,
    pon_path: str | None = None,
    only_contigs: list[str] | None = None,
    threads: int = 8,
    command: str = "",
) -> tuple[Path, Path]:
    """Run single-sample calling; returns (vcf_path, store_path)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = config.reference.fasta
    if not fasta:
        raise ValueError("reference.fasta is not configured")
    if not depth_bedgraph and not cram:
        raise ValueError("Provide a CRAM (self-contained) and/or a depth bedgraph (fast path)")

    all_contigs = read_fai(fasta)
    by_name = {c.name: c for c in all_contigs}
    call_names = [
        n for n in config.reference.call_contigs
        if n in by_name and (only_contigs is None or n in only_contigs)
    ]
    if not call_names:
        raise ValueError("No callable contigs selected")

    # ---- CRAM scan (PE/SR evidence, and depth on the self-contained path)
    insert = InsertModel(mu=0.0, sigma=0.0, n_pairs=0)
    scan_results = None
    if cram:
        from .extract.scan import scan_cram

        insert, scan_results = scan_cram(
            cram, fasta, {n: by_name[n].length for n in call_names},
            config.extract, threads=threads,
        )

    bin_size = config.extract.fast_bin_size if depth_bedgraph else config.extract.bin_size
    bins = GenomeBins([by_name[n] for n in call_names], bin_size)
    tracks = ReferenceTracks(
        build_tracks(fasta, [by_name[n] for n in call_names], bin_size, resources_dir)
    )
    exclude = read_bed(config.reference.exclude_bed) if config.reference.exclude_bed else None

    # ---- depth + masks
    raw, mask = {}, {}
    for name in call_names:
        n_bins = bins.n_bins[name]
        if depth_bedgraph:
            raw[name] = depth_mod.read_depth_bedgraph(depth_bedgraph, name, n_bins, bin_size)
        else:
            raw[name] = scan_results[name].depth
        mask[name] = build_bin_mask(n_bins, bin_size, name, tracks.nfrac(name), exclude)

    autosomes = [n for n in call_names if by_name[n].kind == "autosome"]
    if not autosomes:  # single-contig smoke runs on chrX etc.
        autosomes = call_names
    gc_model = depth_mod.fit_gc_model(
        raw, {n: tracks.gc(n) for n in call_names}, mask, autosomes
    )
    corr = {n: depth_mod.correct_depth(raw[n], tracks.gc(n), gc_model) for n in call_names}
    stats = depth_mod.sample_depth_stats(corr, mask, autosomes)
    mean_depth = stats["mean_depth"]
    if mean_depth <= 0:
        raise ValueError(f"{sample}: zero usable depth")
    ratio = {n: corr[n] / mean_depth for n in call_names}
    # ---- BAF
    baf: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for name in call_names:
        if vaf_bedgraph:
            baf[name] = baf_mod.read_vaf_bedgraph(vaf_bedgraph, name)
        elif snv_vcf:
            baf[name] = baf_mod.read_vaf_vcf(snv_vcf, name, config.extract.baf_min_dp)
        else:
            baf[name] = (
                np.array([], np.uint32), np.array([], np.float32), np.array([], np.uint16)
            )

    ploidy = estimate_sex_ploidy(ratio, mask, by_name)
    karyo = None
    if pon_path and "chrX" in ratio:
        # the X/Y baseline from the measured karyotype (copy number on the
        # panel's sex-informative bins, X heterozygosity, BAF): the threshold
        # inference above reads an XXY as a man (chrX then called as 172
        # fragmented DUPs, a DEL of one X invisible), an XXX as a woman (242)
        karyo = sex_karyotype(sample, bin_size, raw, mask, baf, pon_path)
        sp = sex_ploidy(karyo["karyotype"])
        if sp is not None:
            for name in ratio:
                if by_name[name].kind in ("X", "Y"):
                    ploidy[name] = sp["chrX" if by_name[name].kind == "X" else "chrY"]
        logger.info("%s: karyotype %s%s -> chrX %s, chrY %s", sample, karyo["karyotype"],
                    f' ({karyo["note"]})' if karyo["note"] else "",
                    ploidy.get("chrX"), ploidy.get("chrY"))
    # whole autosomes at 1 or 3 copies (BAF-confirmed): their own baseline
    aneu = autosomal_aneuploidies(
        ratio, mask, baf, bin_size, [n for n in autosomes if by_name[n].kind == "autosome"]
    )
    for name, k in aneu.items():
        ploidy[name] = k
        logger.info("%s: %s at %d copies (whole chromosome) -> baseline %d", sample, name, k, k)

    # ---- wave QC (megabase-scale coverage oscillations surviving GC
    #      correction), on the depth as sequenced: every threshold built on
    #      wave_score (QC alerts, sub-threshold tier, row bounds) is calibrated
    #      on it. With a PoN the corrected depth is scored too (wave_score_pon):
    #      on PMS a founder panel took the median from 0.088 to 0.020.
    from .qc import sample_wave_score

    wave_score, wave_by_contig = sample_wave_score(ratio, mask, bin_size)
    if wave_score is not None:
        logger.info("%s: wave score %.3f", sample, wave_score)
    wave_score_pon = None

    # ---- panel-of-normals normalization (divides out residual systematic bias
    #      and masks population-variable bins)
    if pon_path:
        from .cohort.pon import PanelOfNormals

        pon = PanelOfNormals(pon_path)
        if pon.bin_size != bin_size:
            raise ValueError(
                f"PoN bin size {pon.bin_size} != run bin size {bin_size}"
            )
        n_new_masked = 0
        for name in call_names:
            before = int(mask[name].sum())
            sex, level = pon_profile(by_name[name].kind, ploidy.get(name, 2))
            ratio[name], mask[name] = pon.normalize(
                name, ratio[name], mask[name], sex=sex, level=level
            )
            corr[name] = (ratio[name] * mean_depth).astype(np.float32)
            n_new_masked += int(mask[name].sum()) - before
        logger.info("PoN: applied (chrX %s), %d additional bins masked",
                    pon_profile("X", ploidy.get("chrX", 2))[0], n_new_masked)
        wave_score_pon, _ = sample_wave_score(ratio, mask, bin_size)
        if wave_score_pon is not None:
            logger.info("%s: wave score after PoN %.3f", sample, wave_score_pon)

    # ---- evidence store
    store_path = out_dir / f"{sample}.evidence.zarr"
    writer = EvidenceStoreWriter(store_path, sample, reference_id(fasta), bin_size)
    writer.set_meta(
        mean_depth=mean_depth,
        depth_mad=stats["depth_mad"],
        ploidy={k: int(v) for k, v in ploidy.items()},
        karyotype=karyo["karyotype"] if karyo else None,
        karyotype_note=karyo["note"] if karyo else None,
        aneuploidies={k: int(v) for k, v in aneu.items()} or None,
        depth_source="bedgraph" if depth_bedgraph else "cram",
        insert_mu=insert.mu,
        insert_sigma=insert.sigma,
        insert_n_pairs=insert.n_pairs,
        wave_score=wave_score,
        wave_by_contig={k: round(v, 4) for k, v in wave_by_contig.items()},
        wave_score_pon=wave_score_pon,
        pon=str(pon_path) if pon_path else None,
        baf_source=("vaf-bedgraph" if vaf_bedgraph else
                    f"vcf:{Path(snv_vcf).name} PASS DP>={config.extract.baf_min_dp}"
                    if snv_vcf else None),
    )
    for name in call_names:
        writer.write_depth(name, raw[name], corr[name], mask[name])
        pos, vaf, dp = baf[name]
        if len(pos):
            writer.write_baf(name, pos, vaf, dp)
        if scan_results is not None:
            ev = scan_results[name]
            writer.write_pe(name, **ev.pe)
            writer.write_sr(name, **ev.sr, seqs=ev.sr_seq)
            writer.write_cig(name, **ev.cig)
            if ev.ins:
                writer.write_ins(name, **ev.ins)

    # ---- segmentation, junctions, reconciliation, genotyping
    # sequence coverage estimate (for the junction-support Poisson rate)
    if depth_bedgraph:
        coverage = mean_depth  # bedgraph values are per-base coverage
    else:
        # read-start counts per bin -> coverage = counts * readlen / bin_size
        coverage = mean_depth * (insert.read_length or 150.0) / bin_size
    writer.set_meta(coverage=coverage)
    genotyper = Genotyper(config.genotype, config.segment.max_cn)
    calls: list[CNVCall] = []
    for name in call_names:
        if ploidy[name] == 0:
            continue  # chrY in females
        pos, vaf, _dp = baf[name]
        baf_emit = None
        if len(pos) and config.segment.baf_weight > 0:
            baf_emit = hmm.baf_log_emissions(
                pos, vaf, bins.n_bins[name], bin_size, config.segment.max_cn
            )
        segments, _path = hmm.segment_contig(
            ratio[name], mask[name], ploidy[name],
            config.segment.max_cn, config.segment.transition,
            baf_emit=baf_emit, baf_weight=config.segment.baf_weight,
            bridge_max_gap_bins=config.segment.bridge_max_gap_bins,
            bridge_max_gap_frac=config.segment.bridge_max_gap_frac,
            bridge_max_ratio_diff=config.segment.bridge_max_ratio_diff,
        )
        junctions = []
        if scan_results is not None and insert.n_pairs:
            from .breakpoints.cluster import cluster_cigar_dels
            from .breakpoints.resolve import merge_junction_sets

            ev = scan_results[name]
            junctions = make_junctions(ev.pe, ev.sr, insert, config.breakpoint)
            cig_junctions = cluster_cigar_dels(ev.cig, config.breakpoint)
            junctions = merge_junction_sets(junctions, cig_junctions)
            logger.info(
                "%s %s: %d CIGAR-deletion junctions", sample, name, len(cig_junctions)
            )
            if config.breakpoint.hidden_splits and ev.sr_seq:
                import pysam as _pysam

                from .breakpoints.resolve import recover_junctions

                with _pysam.FastaFile(fasta) as fa:
                    ref_arr = np.frombuffer(
                        fa.fetch(name).upper().encode(), dtype=np.uint8
                    )
                recovered = recover_junctions(ev.sr, ev.sr_seq, ref_arr, config.breakpoint)
                junctions = merge_junction_sets(junctions, recovered)
                logger.info(
                    "%s %s: %d hidden-split junctions recovered (%d DUP)",
                    sample, name, len(recovered),
                    sum(1 for j in recovered if j.svtype == "DUP"),
                )

        contig_calls = reconcile(
            name, segments, junctions, bin_size, ploidy[name], config
        )
        data = ContigData(
            ratio=ratio[name], mask=mask[name],
            baf_pos=pos, baf_vaf=vaf,
            bin_size=bin_size, ploidy=ploidy[name],
            mean_depth=coverage,
            sigma0=hmm.fit_sigma0(ratio[name], mask[name], ploidy[name]),
        )
        for call in contig_calls:
            if call.source != "INSERT":  # INSERT calls are pre-genotyped by the detector
                genotyper.genotype(call, data)

        # sub-insert DELs from the local insert-size shift (no clean junction)
        if scan_results is not None and insert.n_pairs and scan_results[name].ins:
            from .breakpoints.insert_shift import (
                dedupe_against,
                detect_insert_shift_dels,
            )

            ev = scan_results[name]
            ins_calls = detect_insert_shift_dels(
                name, ev.ins["sum"], ev.ins["n"],
                config.extract.insert_bin_size, insert, config.breakpoint,
            )
            ins_calls = dedupe_against(ins_calls, contig_calls)
            logger.info("%s %s: %d insert-shift DEL candidates", sample, name, len(ins_calls))
            contig_calls.extend(ins_calls)
            contig_calls.sort(key=lambda c: (c.start, c.end))
        calls.extend(contig_calls)
        logger.info(
            "%s %s: %d segments, %d junctions -> %d candidates",
            sample, name, len(segments), len(junctions), len(contig_calls),
        )

    vcf_path = out_dir / f"{sample}.vedette.vcf.gz"
    write_vcf(calls, [by_name[n] for n in call_names], sample, vcf_path, command,
              karyotype=karyo["karyotype"] if karyo else None, aneuploidies=aneu)
    logger.info("%s: %d calls -> %s", sample, len(calls), vcf_path)
    return vcf_path, store_path
