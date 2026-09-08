"""Single-sample calling pipeline: evidence -> evidence store + VCF.

v0 implements the *fast path*: depth from a precomputed dense bedgraph and BAF
from a per-site VAF bedgraph (or SNV VCF), CN-HMM segmentation, depth-only
calls. The self-contained CRAM scan (PE/SR + 100 bp depth) plugs into the same
flow via extract.scan and upgrades calls with junction evidence.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from .config import Config
from .extract import baf as baf_mod
from .extract import depth as depth_mod
from .extract.store import EvidenceStoreWriter
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
from .variants import SOURCE_DEPTH, CNVCall
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


def segments_to_calls(
    contig: str,
    segments: list[hmm.Segment],
    bin_size: int,
    ploidy: int,
    config: Config,
) -> list[CNVCall]:
    calls = []
    for seg in segments:
        if seg.cn == ploidy:
            continue
        svtype = "DEL" if seg.cn < ploidy else "DUP"
        call = CNVCall(
            contig=contig,
            start=seg.start_bin * bin_size,
            end=seg.end_bin * bin_size,
            svtype=svtype,
            source=SOURCE_DEPTH,
            cn=seg.cn,
            cn_posterior=seg.posterior,
            mean_ratio=seg.mean_ratio,
            n_bins=seg.n_bins,
            masked_frac=seg.masked_frac,
            cipos=(-bin_size, bin_size),
            ciend=(-bin_size, bin_size),
            precise=False,
        )
        # provisional depth-only genotype (refined by genotype.likelihood once
        # junction/BAF likelihoods are wired in)
        if ploidy == 2:
            if seg.cn == 1 or seg.cn == 3:
                call.gt = (0, 1)
            elif seg.cn == 0 or seg.cn >= 4:
                call.gt = (1, 1)
        else:  # haploid contig
            call.gt = (1, 1)
        call.gq = int(min(99, -10 * np.log10(max(1 - seg.posterior, 1e-10))))

        if seg.n_bins < config.segment.min_bins:
            call.filters.append("NoEvidence")
        if seg.masked_frac > config.segment.max_masked_frac:
            call.filters.append("MaxN")
        if seg.posterior < config.segment.min_cn_posterior:
            call.filters.append("LowPosterior")
        if call.size < config.breakpoint.min_size:
            call.filters.append("MinSize")
        calls.append(call)
    return calls


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
    command: str = "",
) -> tuple[Path, Path]:
    """Run single-sample calling; returns (vcf_path, store_path)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = config.reference.fasta
    if not fasta:
        raise ValueError("reference.fasta is not configured")

    all_contigs = read_fai(fasta)
    by_name = {c.name: c for c in all_contigs}
    call_names = [
        n for n in config.reference.call_contigs
        if n in by_name and (only_contigs is None or n in only_contigs)
    ]
    if not call_names:
        raise ValueError("No callable contigs selected")

    if depth_bedgraph:
        bin_size = config.extract.fast_bin_size
    else:
        raise NotImplementedError(
            "Self-contained CRAM depth extraction lands with extract.scan; "
            "provide --depth-bedgraph for now"
        )

    bins = GenomeBins([by_name[n] for n in call_names], bin_size)
    tracks = ReferenceTracks(
        build_tracks(fasta, [by_name[n] for n in call_names], bin_size, resources_dir)
    )
    exclude = read_bed(config.reference.exclude_bed) if config.reference.exclude_bed else None

    # ---- depth + masks
    raw, mask = {}, {}
    for name in call_names:
        n_bins = bins.n_bins[name]
        raw[name] = depth_mod.read_depth_bedgraph(depth_bedgraph, name, n_bins, bin_size)
        mask[name] = build_bin_mask(
            n_bins, bin_size, name, tracks.nfrac(name), exclude
        )
    autosomes = [n for n in call_names if by_name[n].kind == "autosome"]
    if not autosomes:  # single-contig smoke runs on chrX etc.
        autosomes = call_names
    gc_model = depth_mod.fit_gc_model(
        raw, {n: tracks.gc(n) for n in call_names}, mask, autosomes
    )
    corr = {
        n: depth_mod.correct_depth(raw[n], tracks.gc(n), gc_model) for n in call_names
    }
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
        insert_mu=0.0,
        insert_sigma=0.0,
    )
    for name in call_names:
        writer.write_depth(name, raw[name], corr[name], mask[name])
        pos, vaf, dp = baf[name]
        if len(pos):
            writer.write_baf(name, pos, vaf, dp)

    # ---- segmentation and calls
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
            ratio[name],
            mask[name],
            ploidy[name],
            config.segment.max_cn,
            config.segment.transition,
            baf_emit=baf_emit,
            baf_weight=config.segment.baf_weight,
        )
        contig_calls = segments_to_calls(name, segments, bin_size, ploidy[name], config)
        # annotate BAF summary per call
        for call in contig_calls:
            if len(pos):
                sel = (pos >= call.start) & (pos < call.end)
                het = (vaf[sel] > 0.15) & (vaf[sel] < 0.85)
                call.baf_n = int(het.sum())
                if call.baf_n:
                    call.baf_dev = float(np.mean(np.abs(vaf[sel][het] - 0.5)))
        calls.extend(contig_calls)
        logger.info("%s %s: %d candidate calls", sample, name, len(contig_calls))

    vcf_path = out_dir / f"{sample}.vedette.vcf.gz"
    write_vcf(calls, [by_name[n] for n in call_names], sample, vcf_path, command)
    logger.info("%s: %d calls -> %s", sample, len(calls), vcf_path)
    return vcf_path, store_path
