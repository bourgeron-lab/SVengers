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


def run_call(
    sample: str,
    out_dir: str | Path,
    config: Config,
    resources_dir: str | Path,
    cram: str | None = None,
    depth_bedgraph: str | None = None,
    vaf_bedgraph: str | None = None,
    snv_vcf: str | None = None,
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
    ploidy = estimate_sex_ploidy(ratio, mask, by_name)

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

    # ---- evidence store
    store_path = out_dir / f"{sample}.evidence.zarr"
    writer = EvidenceStoreWriter(store_path, sample, reference_id(fasta), bin_size)
    writer.set_meta(
        mean_depth=mean_depth,
        depth_mad=stats["depth_mad"],
        ploidy={k: int(v) for k, v in ploidy.items()},
        depth_source="bedgraph" if depth_bedgraph else "cram",
        insert_mu=insert.mu,
        insert_sigma=insert.sigma,
        insert_n_pairs=insert.n_pairs,
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
        )
        junctions = []
        if scan_results is not None and insert.n_pairs:
            ev = scan_results[name]
            junctions = make_junctions(ev.pe, ev.sr, insert, config.breakpoint)
            if config.breakpoint.hidden_splits and ev.sr_seq:
                import pysam as _pysam

                from .breakpoints.resolve import merge_junction_sets, recover_junctions

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
            sigma0=hmm.fit_sigma0(ratio[name], mask[name]),
        )
        for call in contig_calls:
            genotyper.genotype(call, data)
        calls.extend(contig_calls)
        logger.info(
            "%s %s: %d segments, %d junctions -> %d candidates",
            sample, name, len(segments), len(junctions), len(contig_calls),
        )

    vcf_path = out_dir / f"{sample}.vedette.vcf.gz"
    write_vcf(calls, [by_name[n] for n in call_names], sample, vcf_path, command)
    logger.info("%s: %d calls -> %s", sample, len(calls), vcf_path)
    return vcf_path, store_path
