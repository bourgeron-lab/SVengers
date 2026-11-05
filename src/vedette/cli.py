"""Vedette command-line interface."""

from __future__ import annotations

import logging
import shlex
import sys
from pathlib import Path
from typing import Optional

import typer

from . import __version__
from .config import load_config

app = typer.Typer(
    name="vedette",
    help="Vedette: multi-evidence germline CNV caller (DEL/DUP) with calibrated filtering.",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"vedette {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None, "--version", callback=_version_callback, is_eager=True
    ),
) -> None:
    pass


@app.command()
def call(
    sample: str = typer.Option(..., help="Sample identifier (matches the CRAM read group)"),
    out_dir: Path = typer.Option(..., help="Output directory (VCF + evidence store)"),
    fasta: Optional[Path] = typer.Option(None, help="Reference FASTA (with .fai)"),
    resources: Path = typer.Option(
        Path("vedette_resources"), help="Writable directory for cached reference tracks"
    ),
    cram: Optional[Path] = typer.Option(None, help="Input CRAM (self-contained path)"),
    depth_bedgraph: Optional[Path] = typer.Option(
        None, help="Precomputed dense depth bedgraph (fast path)"
    ),
    vaf_bedgraph: Optional[Path] = typer.Option(
        None, help="Precomputed per-site VAF bedgraph (fast path BAF)"
    ),
    snv_vcf: Optional[Path] = typer.Option(None, help="Single-sample SNV VCF for BAF"),
    pon: Optional[Path] = typer.Option(
        None, help="Panel-of-normals (.zarr) for residual depth normalization"
    ),
    exclude_bed: Optional[Path] = typer.Option(None, help="Exclusion regions BED"),
    contig: Optional[list[str]] = typer.Option(
        None, "--contig", help="Restrict to contig(s), repeatable (smoke tests)"
    ),
    config_file: Optional[Path] = typer.Option(None, "--config", help="TOML configuration"),
) -> None:
    """Single-sample CNV calling: evidence extraction + segmentation -> VCF + evidence store."""
    from .call import run_call

    overrides: dict = {"reference": {}}
    if fasta:
        overrides["reference"]["fasta"] = str(fasta)
    if exclude_bed:
        overrides["reference"]["exclude_bed"] = str(exclude_bed)
    config = load_config(config_file, **overrides)

    vcf, store = run_call(
        sample=sample,
        out_dir=out_dir,
        config=config,
        resources_dir=resources,
        cram=str(cram) if cram else None,
        depth_bedgraph=str(depth_bedgraph) if depth_bedgraph else None,
        vaf_bedgraph=str(vaf_bedgraph) if vaf_bedgraph else None,
        snv_vcf=str(snv_vcf) if snv_vcf else None,
        pon_path=str(pon) if pon else None,
        only_contigs=contig if contig else None,
        command=shlex.join(sys.argv),
    )
    typer.echo(f"VCF:   {vcf}")
    typer.echo(f"Store: {store}")


@app.command()
def ped(
    pedigree: Path = typer.Argument(..., help="6-column pedigree TSV (header optional)"),
) -> None:
    """Parse and summarize a pedigree file (Mendelian annotation arrives with the cohort stage)."""
    from .ped import parse_ped

    p = parse_ped(pedigree)
    typer.echo(
        f"{len(p.individuals)} individuals, {len(p.families)} families, "
        f"{len(p.trios())} trios, {len(p.founders())} founders"
    )


manifest_app = typer.Typer(help="Cohort sample manifests")
app.add_typer(manifest_app, name="manifest")


@manifest_app.command("discover")
def manifest_discover(
    root: Path = typer.Argument(..., help="Samples root directory"),
    out: Path = typer.Option(..., help="Output manifest TSV"),
    layout: str = typer.Option("ghfc", help="Directory layout (currently: ghfc)"),
    pedigree: Optional[Path] = typer.Option(
        None, help="Restrict to samples of this pedigree"
    ),
) -> None:
    """Build a manifest by walking a known samples directory layout."""
    from .manifest import discover_ghfc
    from .ped import parse_ped

    if layout != "ghfc":
        raise typer.BadParameter(f"Unknown layout {layout!r}")
    ids = list(parse_ped(pedigree).individuals) if pedigree else None
    manifest = discover_ghfc(root, sample_ids=ids)
    manifest.write(out)
    typer.echo(f"{len(manifest.samples)} samples -> {out}")


@manifest_app.command("validate")
def manifest_validate(manifest: Path = typer.Argument(...)) -> None:
    """Check that every CRAM and index in a manifest exists."""
    from .manifest import Manifest

    problems = Manifest.read(manifest).validate()
    for problem in problems:
        typer.echo(problem, err=True)
    if problems:
        raise typer.Exit(code=1)
    typer.echo("OK")


def _not_yet(name: str) -> None:
    typer.echo(f"`vedette {name}` is not implemented yet in this development snapshot.", err=True)
    raise typer.Exit(code=2)


@app.command()
def filter(
    vcf: Path = typer.Argument(..., help="Single-sample Vedette VCF"),
    store: Path = typer.Option(..., help="Matching evidence store (.zarr)"),
    out: Path = typer.Option(..., help="Output filtered VCF (.vcf.gz)"),
    segdup_bed: Optional[Path] = typer.Option(None, help="Segmental duplications BED"),
    lcr_bed: Optional[Path] = typer.Option(None, help="LCR/self-chain BED"),
    exclude_bed: Optional[Path] = typer.Option(None, help="Exclusion regions BED"),
    popcov_bw: Optional[Path] = typer.Option(None, help="Population coverage stdev bigWig"),
    model_dir: Optional[Path] = typer.Option(
        None, help="Trained filter model directory (enables SQ scoring)"
    ),
    features_out: Optional[Path] = typer.Option(
        None, help="Also export the per-call feature matrix (TSV)"
    ),
    artifact_bed: Optional[Path] = typer.Option(
        None, help="Recurrent-artifact BED (from `vedette pon artifacts`)"
    ),
    config_file: Optional[Path] = typer.Option(None, "--config", help="TOML configuration"),
) -> None:
    """Tiered false-positive filtering (hard filters + calibrated SQ score)."""
    from .filter.pipeline import run_filter

    overrides: dict = {"reference": {}}
    if segdup_bed:
        overrides["reference"]["segdup_bed"] = str(segdup_bed)
    if lcr_bed:
        overrides["reference"]["lcr_bed"] = str(lcr_bed)
    if exclude_bed:
        overrides["reference"]["exclude_bed"] = str(exclude_bed)
    if popcov_bw:
        overrides["reference"]["pop_cov_stdev_bw"] = str(popcov_bw)
    config = load_config(config_file, **overrides)
    counts = run_filter(
        vcf, store, out, config,
        model_dir=model_dir, features_out=features_out, artifact_bed=artifact_bed,
    )
    typer.echo(
        f"{counts['total']} calls: {counts['pass']} PASS, "
        f"{counts['hard_filtered']} hard-filtered, {counts['low_sq']} LowSQ -> {out}"
    )


@app.command()
def tensors(
    vcf: Path = typer.Argument(..., help="Filtered single-sample VCF (from `vedette filter`)"),
    store: Path = typer.Option(..., help="Matching evidence store (.zarr)"),
    out: Path = typer.Option(..., help="Output tensors (.npz)"),
    min_size: int = typer.Option(1000, help="Smallest call kept (bp)"),
    all_filters: bool = typer.Option(
        False, help="Keep hard-filtered calls too (default: PASS and LowSQ only)"
    ),
) -> None:
    """Per-call multi-channel tensors for the CNN second opinion, from the store."""
    from .filter.tensor import CHANNELS, L, run_tensors

    n = run_tensors(vcf, store, out, min_size=min_size,
                    filters=None if all_filters else ("PASS", "LowSQ"))
    typer.echo(f"{n} calls x {len(CHANNELS)} channels x {L} points -> {out}")


@app.command()
def merge(
    vcfs: list[Path] = typer.Argument(..., help="Filtered single-sample VCFs (or one @filelist)"),
    out: Path = typer.Option(..., help="Output sites-only VCF (.vcf.gz)"),
    fasta: Path = typer.Option(..., help="Reference FASTA (for the contig header)"),
    min_gq: int = typer.Option(0, help="Minimum GQ for a call to nominate a site"),
    include_filtered: bool = typer.Option(
        False, help="Also cluster non-PASS calls (default: PASS/LowPosterior only)"
    ),
    rescue_min_sq: float = typer.Option(
        1.0, help="LowSQ calls at or above this SQ join clusters as sub-threshold members (<0 disables)"
    ),
    rescue_min_samples: int = typer.Option(
        2, help="Distinct sub-threshold carriers required to rescue a site with no PASS member"
    ),
    rescue_total_sq: float = typer.Option(
        8.0, help="Summed best-SQ across carriers required to rescue such a site"
    ),
    jobs: int = typer.Option(0, help="Worker processes (0 = every allowed CPU)"),
) -> None:
    """Cross-sample clustering of filtered calls into cohort sites."""
    from .cohort.merge import run_merge
    from .parallel import resolve_jobs
    from .reference import read_fai

    paths: list[Path] = []
    for v in vcfs:
        if str(v).startswith("@"):
            paths.extend(Path(p) for p in Path(str(v)[1:]).read_text().split())
        else:
            paths.append(v)
    n = run_merge(
        paths, read_fai(fasta), out, min_gq=min_gq, pass_only=not include_filtered,
        rescue_min_sq=None if rescue_min_sq < 0 else rescue_min_sq,
        rescue_min_samples=rescue_min_samples, rescue_total_sq=rescue_total_sq,
        jobs=resolve_jobs(jobs),
    )
    typer.echo(f"{len(paths)} VCFs -> {n} sites -> {out}")


@app.command()
def genotype(
    sites: Path = typer.Argument(..., help="Sites VCF from `vedette merge`"),
    store: Path = typer.Option(..., help="This sample's evidence store (.zarr)"),
    out: Path = typer.Option(..., help="Output per-sample genotype table (TSV)"),
    config_file: Optional[Path] = typer.Option(None, "--config", help="TOML configuration"),
) -> None:
    """Force-genotype cohort sites in ONE sample from its evidence store (no CRAM)."""
    from .cohort.genotyper import run_genotype

    n = run_genotype(sites, store, out, load_config(config_file))
    typer.echo(f"{n} site genotypes -> {out}")


@app.command()
def mosaic(
    store: Optional[Path] = typer.Argument(None, help="One sample's evidence store (.zarr)"),
    out: Optional[Path] = typer.Option(None, help="Output TSV (with STORE)"),
    run_dir: Optional[Path] = typer.Option(
        None, "--run-dir",
        help="Scan every samples/<S>/<S>.evidence.zarr of a run, writing <S>.mosaic.tsv beside it",
    ),
) -> None:
    """Chromosomal-scale mosaic CNVs (fractional copy number confirmed by BAF), from the store."""
    from .extract.store import EvidenceStore
    from .mosaic import scan_store, write_tsv

    if run_dir is not None:
        n_samples = n_calls = n_carriers = 0
        for d in sorted((run_dir / "samples").iterdir()):
            path = d / f"{d.name}.evidence.zarr"
            if not path.exists():
                continue
            try:
                calls = scan_store(EvidenceStore(path))
            except Exception as exc:  # one unreadable store must not stop the batch
                typer.echo(f"{d.name}: skipped ({exc})", err=True)
                continue
            write_tsv(calls, d / f"{d.name}.mosaic.tsv")
            n_samples += 1
            n_calls += len(calls)
            n_carriers += any(c.mosaic for c in calls)
        typer.echo(f"{n_samples} samples scanned, {n_calls} segments, "
                   f"{n_carriers} samples with a mosaic segment")
        return
    if store is None or out is None:
        raise typer.BadParameter("give STORE and --out, or --run-dir")
    calls = scan_store(EvidenceStore(store))
    write_tsv(calls, out)
    typer.echo(f"{len(calls)} segments ({sum(c.mosaic for c in calls)} mosaic) -> {out}")


@app.command()
def wave(
    stores: list[Path] = typer.Argument(None, help="Evidence stores (.zarr)"),
    run_dir: Optional[Path] = typer.Option(
        None, "--run-dir", help="Every samples/<S>/<S>.evidence.zarr of a run"
    ),
    force: bool = typer.Option(False, "--force", help="Recompute even where a score exists"),
    jobs: int = typer.Option(1, "--jobs", help="Stores processed in parallel"),
) -> None:
    """Write the wave QC score into stores that lack it (called before schema 0.4)."""
    from concurrent.futures import ProcessPoolExecutor

    from .qc import WAVE_WARN, backfill_wave

    paths = list(stores or [])
    if run_dir is not None:
        paths += [d / f"{d.name}.evidence.zarr" for d in sorted((run_dir / "samples").iterdir())
                  if (d / f"{d.name}.evidence.zarr").exists()]
    if not paths:
        raise typer.BadParameter("give STORE(s) or --run-dir")
    n_written = n_wavy = n_failed = 0
    with ProcessPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {pool.submit(backfill_wave, p, force): p for p in paths}
        for fut, p in futures.items():
            try:
                score, written = fut.result()
            except Exception as exc:  # one unreadable store must not stop the batch
                typer.echo(f"{p}: skipped ({exc})", err=True)
                n_failed += 1
                continue
            n_written += written
            n_wavy += score is not None and score >= WAVE_WARN
    typer.echo(f"{len(paths)} stores: {n_written} scores written, {n_wavy} >= {WAVE_WARN} "
               f"(wavy), {n_failed} failed")


@app.command()
def baf(
    store: Optional[Path] = typer.Argument(None, help="Evidence store (.zarr)"),
    vcf: Optional[Path] = typer.Option(None, "--vcf", help="Single-sample SNV VCF of the store's sample"),
    manifest: Optional[Path] = typer.Option(
        None, "--manifest", help="TSV of store<TAB>vcf pairs, one per line (instead of STORE --vcf)"
    ),
    min_dp: int = typer.Option(10, "--min-dp", help="Minimum site depth (FORMAT/DP, else the sum of AD)"),
    force: bool = typer.Option(False, "--force", help="Rewrite stores whose BAF already comes from a VCF"),
    jobs: int = typer.Option(1, "--jobs", help="Stores processed in parallel"),
) -> None:
    """Rewrite the BAF of evidence stores from SNV VCFs (PASS, no QUAL threshold), without the CRAM."""
    from concurrent.futures import ProcessPoolExecutor

    from .extract.baf import backfill_baf

    pairs: list[tuple[Path, Path]] = []
    if store is not None:
        if vcf is None:
            raise typer.BadParameter("STORE needs --vcf")
        pairs.append((store, vcf))
    if manifest is not None:
        for line in manifest.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                s, v = line.split("\t")[:2]
                pairs.append((Path(s), Path(v)))
    if not pairs:
        raise typer.BadParameter("give STORE --vcf VCF or --manifest")
    n_written = n_failed = 0
    with ProcessPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {pool.submit(backfill_baf, s, str(v), min_dp, force): s for s, v in pairs}
        for fut, s in futures.items():
            try:
                n_sites, n_het, written = fut.result()
            except Exception as exc:  # one unreadable store or VCF must not stop the batch
                typer.echo(f"{s}: skipped ({exc})", err=True)
                n_failed += 1
                continue
            n_written += written
            if written:
                typer.echo(f"{s}: {n_sites} sites, {n_het} heterozygous")
    typer.echo(f"{len(pairs)} stores: {n_written} rewritten, {n_failed} failed")


@app.command()
def regenotype(
    vcf: Path = typer.Argument(..., help="The sample's raw calls VCF (<S>.vedette.vcf.gz)"),
    store: Path = typer.Option(..., help="The sample's evidence store"),
    out: Path = typer.Option(..., help="Output VCF (bgzipped, indexed)"),
    config_file: Optional[Path] = typer.Option(None, "--config", help="TOML configuration"),
) -> None:
    """Recompute GT, GQ, CN, CNQ, RD, BAFN, BAFDEV of a sample's calls from its store (after `vedette baf`)."""
    from .genotype.regenotype import regenotype_calls

    c = regenotype_calls(vcf, store, out, load_config(config_file))
    typer.echo(f"{c['calls']} calls: {c['regenotyped']} re-genotyped ({c['gt_changed']} GT changed), "
               f"{c['unchanged']} unchanged (INSERT or no depth) -> {out}")


@app.command("site-map")
def site_map(
    old: Path = typer.Argument(..., help="Cohort VCF of the previous run"),
    new: Path = typer.Argument(..., help="Cohort VCF of the new run"),
    out: Path = typer.Option(..., help="Output TSV: old site -> new site(s), with the relation"),
    min_ro: float = typer.Option(0.5, help="Reciprocal overlap for two sites to be the same event"),
) -> None:
    """Correspondence of cohort sites between two runs (same, moved, split, merged, lost, new)."""
    from .cohort.sitemap import map_sites, write_map

    rows, counts = map_sites(old, new, min_ro)
    write_map(rows, out)
    order = ("same", "moved", "split", "merged", "complex", "lost", "new")
    typer.echo(", ".join(f"{counts.get(k, 0)} {k}" for k in order) + f" -> {out}")


@app.command("cohort-filter")
def cohort_filter(
    sites: Path = typer.Argument(..., help="Sites VCF from `vedette merge`"),
    genotypes: list[Path] = typer.Option(
        ..., "--genotypes", help="Per-sample genotype TSVs (repeatable, or one @filelist)"
    ),
    out: Path = typer.Option(..., help="Output multi-sample cohort VCF (.vcf.gz)"),
    pedigree: Optional[Path] = typer.Option(
        None, help="Pedigree TSV (founder-based HWE)"
    ),
    artifact_af: float = typer.Option(0.02, help="AF above which HWE failures are flagged"),
    carrier_frac_max: float = typer.Option(
        0.85, help="HighCarrierFrac: carrier fraction at or above which a site is flagged (needs >=4 genotyped samples)"
    ),
    cn_discordant_frac: float = typer.Option(
        0.5, help="CNDiscordant: fraction of carriers whose re-genotyped CN contradicts the SV type"
    ),
    known_sites: Optional[Path] = typer.Option(
        None, help="Reference cohort VCF (e.g. a large `cohort-filter` output): annotates REF_AF/REF_CARRIER_FRAC and flags CommonRef"
    ),
    known_common_frac: float = typer.Option(
        0.5, help="CommonRef: reference-cohort carrier fraction at or above which the site is flagged"
    ),
    known_artifacts: bool = typer.Option(
        True, "--known-artifacts/--no-known-artifacts",
        help="With --known-sites: a site common here (AF > --artifact-af) whose reference site is CohortArtifact is CohortArtifact too (INFO REF_ARTIFACT)",
    ),
    jobs: int = typer.Option(0, help="Worker processes (0 = every allowed CPU)"),
) -> None:
    """Assemble genotype tables into a cohort VCF with AF/HWE/artifact metrics."""
    from .cohort.metrics import run_cohort_filter
    from .parallel import resolve_jobs
    from .ped import parse_ped

    paths: list[Path] = []
    for g in genotypes:
        if str(g).startswith("@"):
            paths.extend(Path(p) for p in Path(str(g)[1:]).read_text().split())
        else:
            paths.append(g)
    counts = run_cohort_filter(
        sites, paths, out,
        pedigree=parse_ped(pedigree) if pedigree else None,
        artifact_af=artifact_af,
        carrier_frac_max=carrier_frac_max,
        cn_discordant_frac=cn_discordant_frac,
        known_sites=known_sites,
        known_common_frac=known_common_frac,
        known_artifacts=known_artifacts,
        jobs=resolve_jobs(jobs),
    )
    typer.echo(
        f"{counts['sites']} sites, {counts['pass']} PASS, {counts['artifact']} CohortArtifact "
        f"({counts['ref_artifact']} from the reference), "
        f"{counts['high_carrier']} HighCarrierFrac, {counts['cn_discordant']} CNDiscordant, "
        f"{counts['common_ref']} CommonRef, {counts['chry_no_depth']} ChrYNoDepth -> {out}"
    )


@app.command("ped-annotate")
def ped_annotate(
    cohort_vcf: Path = typer.Argument(..., help="Multi-sample cohort VCF"),
    pedigree: Path = typer.Option(..., help="Pedigree TSV"),
    out_report: Path = typer.Option(..., help="Per-trio-site report TSV"),
    out_summary: Optional[Path] = typer.Option(None, help="Per-trio summary TSV"),
    jobs: int = typer.Option(0, help="Worker processes (0 = every allowed CPU)"),
) -> None:
    """Mendelian consistency, inheritance and de novo annotation on a cohort VCF."""
    from .cohort.mendel import run_ped_annotate
    from .parallel import resolve_jobs
    from .ped import parse_ped

    stats = run_ped_annotate(
        cohort_vcf, parse_ped(pedigree), out_report, out_summary,
        jobs=resolve_jobs(jobs),
    )
    for st in stats:
        tr = f"{st.transmission_rate:.3f}" if st.n_transmit_opportunity else "NA"
        typer.echo(
            f"{st.family_id}/{st.child}: {st.n_genotyped} genotyped, "
            f"{st.n_child_carrier} carrier, {st.n_denovo} de novo, "
            f"MVR {st.violation_rate:.3f}, transmission {tr}"
        )


pon_app = typer.Typer(help="Panel of normals: depth profile and artifact catalogue")
app.add_typer(pon_app, name="pon")


@pon_app.command("build")
def pon_build(
    stores: list[Path] = typer.Argument(
        ..., help="Evidence stores (.zarr), or one @filelist of paths"
    ),
    out: Path = typer.Option(..., help="Output PoN (.zarr)"),
    min_samples: int = typer.Option(8, help="Minimum samples per profile"),
    denoise_k: int = typer.Option(
        0, help="Also fit k SVD components for residual denoising (0 = off; 5-20 typical)"
    ),
) -> None:
    """Build the per-bin median+MAD depth profile from evidence stores (sex-stratified X/Y)."""
    from .cohort.pon import build_pon

    paths: list[Path] = []
    for s in stores:
        if str(s).startswith("@"):
            paths.extend(Path(p) for p in Path(str(s)[1:]).read_text().split())
        else:
            paths.append(s)
    build_pon(paths, out, min_samples=min_samples, denoise_k=denoise_k)
    typer.echo(f"PoN ({len(paths)} samples) -> {out}")


@pon_app.command("artifacts")
def pon_artifacts(
    cohort_vcf: Path = typer.Argument(..., help="Cohort VCF from `vedette cohort-filter`"),
    out: Path = typer.Option(..., help="Output artifact BED"),
    min_af: float = typer.Option(0.03, help="Founder AF above which HWE failures are artifacts"),
    max_hwe_p: float = typer.Option(1e-6),
    min_gtr_disc: float = typer.Option(0.5),
) -> None:
    """Derive the recurrent-artifact BED from a cohort VCF."""
    from .cohort.pon import artifacts_bed

    n = artifacts_bed(
        cohort_vcf, out, min_af=min_af, max_hwe_p=max_hwe_p, min_gtr_disc=min_gtr_disc
    )
    typer.echo(f"{n} artifact loci -> {out}")


bench_app = typer.Typer(help="Benchmark harness (truth-set evaluation, caller comparison, trios)")
app.add_typer(bench_app, name="bench")


def _stratifier(segdup_bed, lcr_bed, tr_bed):
    from .bench.truth import Stratifier

    return Stratifier(
        segdup_bed=str(segdup_bed) if segdup_bed else "",
        lcr_bed=str(lcr_bed) if lcr_bed else "",
        tr_bed=str(tr_bed) if tr_bed else "",
    )


@bench_app.command("truth")
def bench_truth(
    calls: Path = typer.Argument(..., help="Calls VCF to evaluate"),
    truth: Path = typer.Option(..., help="Truth SV VCF (e.g. GIAB v5.0q stvar)"),
    bed: Path = typer.Option(..., help="Confident-regions BED"),
    out_dir: Path = typer.Option(..., help="Output directory"),
    fasta: Optional[Path] = typer.Option(None, help="Reference FASTA (optional)"),
    segdup_bed: Optional[Path] = typer.Option(None),
    lcr_bed: Optional[Path] = typer.Option(None),
    tr_bed: Optional[Path] = typer.Option(None),
    plot: bool = typer.Option(False, help="Also write F1-by-size figure (PNG)"),
) -> None:
    """Stratified Truvari evaluation of one callset against a truth set."""
    from .bench.truth import evaluate

    rows = evaluate(
        calls, truth, bed, out_dir, reference=fasta,
        stratifier=_stratifier(segdup_bed, lcr_bed, tr_bed),
    )
    for r in rows:
        if r["stratum_kind"] == "all":
            typer.echo(
                f"{r['svtype']}: precision {r['precision']} recall {r['recall']} f1 {r['f1']}"
            )
    if plot:
        from .bench.plots import plot_size_stratified

        plot_size_stratified(rows, Path(out_dir) / "f1_by_size.png")
    typer.echo(f"table -> {out_dir}/eval_stratified.tsv")


@bench_app.command("compare")
def bench_compare(
    entry: list[str] = typer.Option(
        ..., "--entry",
        help="Caller entry name=kind=path (kind: vcf | vcf-all | wisecondorx), repeatable",
    ),
    truth: Path = typer.Option(...),
    bed: Path = typer.Option(...),
    fasta: Path = typer.Option(..., help="Reference FASTA (contig header + truvari)"),
    out_dir: Path = typer.Option(...),
    segdup_bed: Optional[Path] = typer.Option(None),
    lcr_bed: Optional[Path] = typer.Option(None),
    tr_bed: Optional[Path] = typer.Option(None),
    plot: bool = typer.Option(False, help="Also write comparison figure (PNG)"),
) -> None:
    """Run the same truth evaluation over several callers; combined table."""
    from .bench.compare import run_compare
    from .reference import read_fai

    entries = []
    for e in entry:
        parts = e.split("=", 2)
        if len(parts) != 3 or parts[1] not in ("vcf", "vcf-all", "wisecondorx"):
            raise typer.BadParameter(
                f"--entry must be name=vcf|vcf-all|wisecondorx=path, got {e!r}"
            )
        entries.append(tuple(parts))
    rows = run_compare(
        entries, truth, bed, read_fai(fasta), out_dir,
        reference=fasta, stratifier=_stratifier(segdup_bed, lcr_bed, tr_bed),
    )
    for r in rows:
        if r["stratum_kind"] == "all":
            typer.echo(
                f"{r['caller']:>14} {r['svtype']}: precision {r['precision']} "
                f"recall {r['recall']} f1 {r['f1']}"
            )
    if plot and rows:
        from .bench.plots import plot_size_stratified

        plot_size_stratified(rows, Path(out_dir) / "comparison_f1_by_size.png")


@bench_app.command("dup-truth")
def bench_dup_truth(
    truth: Path = typer.Argument(..., help="INS-represented truth VCF (e.g. GIAB v5.0q stvar)"),
    fasta: Path = typer.Option(..., help="Reference FASTA"),
    out: Path = typer.Option(..., help="Output tandem-duplication truth subset (.vcf.gz)"),
    min_identity: float = typer.Option(0.80, help="Flank identity to call an INS tandem"),
) -> None:
    """Isolate the tandem-duplication subset of an INS truth (fair DUP denominator)."""
    from .bench.dupspace import make_tandem_truth

    counts = make_tandem_truth(truth, fasta, out, min_identity=min_identity)
    typer.echo(
        f"{counts['tandem']} tandem / {counts['ins_total']} INS -> {out} "
        f"({counts['other_ins']} insertions non-tandem écartées)"
    )


@bench_app.command("trios")
def bench_trios(
    cohort_vcf: Path = typer.Argument(..., help="Multi-sample cohort VCF"),
    pedigree: Path = typer.Option(...),
    out: Path = typer.Option(..., help="Output sweep TSV"),
    field: str = typer.Option("GQ", help="FORMAT field to threshold (GQ or SQ)"),
    plot: bool = typer.Option(False, help="Also write MVR/transmission figure (PNG)"),
    pass_only: bool = typer.Option(False, help="Restrict to FILTER=PASS sites"),
    exclude_rescue: bool = typer.Option(
        False, help="Also drop INFO/RESCUE sites (recurrent sub-threshold)"
    ),
    jobs: int = typer.Option(0, help="Worker processes (0 = every allowed CPU)"),
) -> None:
    """MVR and transmission rate as a function of the quality threshold."""
    from .bench.trios import sweep
    from .parallel import resolve_jobs
    from .ped import parse_ped

    rows = sweep(
        cohort_vcf, parse_ped(pedigree), field=field, out_tsv=out,
        pass_only=pass_only, exclude_rescue=exclude_rescue,
        jobs=resolve_jobs(jobs),
    )
    if plot:
        from .bench.plots import plot_trio_sweep

        plot_trio_sweep(rows, Path(out).with_suffix(".png"))
    typer.echo(f"sweep -> {out}")


@app.command()
def report(
    run_dir: Path = typer.Argument(..., help="Cohort run directory (Nextflow wrapper layout)"),
    pedigree: Optional[Path] = typer.Option(None, help="Pedigree TSV (enables family sections)"),
    out: Optional[Path] = typer.Option(None, help="Output HTML (default: <run_dir>/report.html)"),
    cohort_name: str = typer.Option("", help="Display name (default: run dir name)"),
    no_family_workspaces: bool = typer.Option(
        False, help="Skip writing families/<FID>/ VCF subsets"
    ),
    min_denovo_gq: int = typer.Option(20, help="GQ threshold for reported de novo candidates"),
    genes_gff: Optional[Path] = typer.Option(
        None, help="Tabix-indexed Ensembl GFF: adds a gene column to site tables"
    ),
    ndd_tsv: Optional[Path] = typer.Option(
        None, help="GeneTrek TSV: NDD tiers feed the impact score of site tables"
    ),
    clingen_genes: Optional[Path] = typer.Option(
        None, help="ClinGen_gene_curation_list_GRCh38.tsv: gene HI/TS scores for the impact score"
    ),
    clingen_regions: Optional[Path] = typer.Option(
        None, help="ClinGen_region_curation_list_GRCh38.tsv: recurrent and other dosage regions (15q11.2q13, 16p11.2, 22q11.2...) for the impact score"
    ),
    segdup_bed: Optional[Path] = typer.Option(
        None, help="Segmental-duplication BED: adds a segdup overlap column to the family events pages"
    ),
    lcr_bed: Optional[Path] = typer.Option(
        None, help="LCR/self-chain BED: adds an LCR overlap column to the family events pages"
    ),
    max_denovo_per_child: int = typer.Option(
        20, help="Children above this many de novo candidates are flagged as technical outliers and left out of family tables"
    ),
    family_events: bool = typer.Option(
        False, "--family-events",
        help="Also write families/<FID>/events.html: every DEL/DUP of the family in one page with chromosome/type filters and inheritance (works with or without family workspaces)",
    ),
    min_carrier_gq: int = typer.Option(
        20, help="Minimum GQ for a genotype to count as a carrier on the events pages (re-genotyping emits low-GQ carriers readily: apparent de novo rate 14-17% below GQ 20 vs 1-2% above)",
    ),
    mosaic_dir: Optional[Path] = typer.Option(
        None, "--mosaic-dir",
        help="Directory of <S>.mosaic.tsv (vedette mosaic); default samples/<S>/ of the run",
    ),
    no_collapse_events: bool = typer.Option(
        False, "--no-collapse-events",
        help="Keep every cohort site as its own row on the events pages (default: merge overlapping same-type events carried by exactly the same members)",
    ),
    events_pass_only: bool = typer.Option(
        False, "--events-pass-only",
        help="Restrict the family events pages to PASS sites (cohort scale: non-PASS artifacts are carried by nearly every family and make the pages unusable)",
    ),
    export_nonpass_max_carrier_frac: float = typer.Option(
        0.20, "--export-nonpass-max-carrier-frac",
        help="With --events-pass-only: the family export (<FID>.vedette.tsv) still holds the non-PASS sites carried by at most this fraction of the cohort (0 = none)",
    ),
    low_quality_prefix: list[str] = typer.Option(
        [], "--low-quality-prefix",
        help="Sample-ID prefix of older, lower-quality genomes (repeatable; IP at the GHFC): no sub-threshold rows for them",
    ),
    jobs: int = typer.Option(0, help="Worker processes (0 = every allowed CPU)"),
    allow_missing_stores: bool = typer.Option(
        False, "--allow-missing-stores",
        help="Report even when no sample has an evidence store under samples/<S>/ (no karyotype section, wave QC or depth-based bounds)",
    ),
) -> None:
    """Self-contained HTML cohort report + per-family workspaces."""
    from .parallel import resolve_jobs
    from .ped import parse_ped
    from .report import generate_report

    try:
        path = generate_report(
            run_dir,
            parse_ped(pedigree) if pedigree else None,
            out_html=out,
            cohort_name=cohort_name,
            make_family_workspaces=not no_family_workspaces,
            min_denovo_gq=min_denovo_gq,
            genes_gff=genes_gff,
            ndd_tsv=ndd_tsv,
            clingen_genes=clingen_genes,
            clingen_regions=clingen_regions,
            segdup_bed=segdup_bed,
            lcr_bed=lcr_bed,
            max_denovo_per_child=max_denovo_per_child,
            family_events=family_events,
            events_pass_only=events_pass_only,
            export_nonpass_max_cf=export_nonpass_max_carrier_frac,
            collapse_family_events=not no_collapse_events,
            min_carrier_gq=min_carrier_gq,
            mosaic_dir=mosaic_dir,
            low_quality_prefixes=tuple(low_quality_prefix),
            jobs=resolve_jobs(jobs),
            require_stores=not allow_missing_stores,
        )
    except FileNotFoundError as e:  # no cohort VCF or no evidence store: no traceback
        typer.secho(f"erreur : {e}", err=True, fg="red")
        raise typer.Exit(1) from None
    typer.echo(f"rapport -> {path}")


@app.command()
def label(
    features: Path = typer.Argument(..., help="Feature TSV from `vedette filter --features-out`"),
    calls_vcf: Path = typer.Option(..., help="The (unfiltered) calls VCF the features came from"),
    truth_vcf: Path = typer.Option(..., help="Truth SV VCF (e.g. GIAB v5.0q stvar)"),
    truth_bed: Path = typer.Option(..., help="Confident-regions BED"),
    fasta: Path = typer.Option(..., help="Reference FASTA"),
    out: Path = typer.Option(..., help="Output labeled feature TSV"),
    truvari_dir: Path = typer.Option(
        Path("truvari_bench"), help="Truvari output directory (recreated)"
    ),
    parent_mode: bool = typer.Option(
        False, help="Transmission-based semi-supervision (sample without its own truth)"
    ),
    recurrent_fp_vcf: Optional[Path] = typer.Option(
        None, help="FP calls of the truth-bearing sample (required with --parent-mode)"
    ),
) -> None:
    """Label calls TP/FP against a truth set (Truvari bench) for filter training."""
    from .train.label import label_features, run_truvari_bench

    run_truvari_bench(calls_vcf, truth_vcf, truth_bed, fasta, truvari_dir)
    counts = label_features(
        features, truvari_dir, out,
        parent_mode=parent_mode, recurrent_fp_vcf=recurrent_fp_vcf,
    )
    typer.echo(
        f"{counts['tp']} TP, {counts['fp']} FP, {counts.get('ambiguous', 0)} ambiguous, "
        f"{counts['unlabeled']} unlabeled -> {out}"
    )


@app.command()
def labels(
    from_genetics_viz: Path = typer.Option(
        ..., "--from-genetics-viz", help="genetics-viz validations/vedette.tsv (curators' validations)"
    ),
    out: Path = typer.Option(..., help="Output TSV: one current validation per (family, variant, sample)"),
    model_version: Optional[str] = typer.Option(
        None, "--model-version", help="Keep only validations of calls made with this model (e.g. v5n)"
    ),
    features: Optional[Path] = typer.Option(
        None, help="Feature TSV of one sample (`vedette filter --features-out`) to label for training"
    ),
    sample: Optional[str] = typer.Option(None, help="With --features: the sample those features belong to"),
    labeled_out: Optional[Path] = typer.Option(None, help="With --features: labeled feature TSV"),
) -> None:
    """Curators' genetics-viz validations as labels: present = 1, absent = 0 (uncertain, different: none)."""
    import statistics

    from .train.validations import (
        breakpoint_offsets,
        label_features,
        read_genetics_viz,
        write_labels,
    )

    vals, counts = read_genetics_viz(from_genetics_viz)
    if model_version:
        vals = [v for v in vals if v.model_version == model_version]
    write_labels(vals, out)
    by = {k: sum(v.validation == k for v in vals) for k in ("present", "absent", "uncertain", "different")}
    typer.echo(f"{counts['rows']} rows ({counts['ignored']} ignored, {counts['superseded']} superseded, "
               f"{counts['unreadable']} unreadable) -> {len(vals)} validations: "
               + ", ".join(f"{n} {k}" for k, n in by.items()) + f" -> {out}")
    offs = breakpoint_offsets(vals)
    for side in ("start", "end"):
        if offs[side]:
            a = [abs(x) for x in offs[side]]
            typer.echo(f"curated {side}: n={len(a)}, median |offset| {statistics.median(a):.0f} bp, "
                       f"<= 1 kb: {sum(x <= 1000 for x in a)}")
    if features is not None:
        if not sample or labeled_out is None:
            raise typer.BadParameter("--features needs --sample and --labeled-out")
        c = label_features(features, vals, sample, labeled_out)
        typer.echo(f"{c['tp']} TP, {c['fp']} FP, {c['unlabeled']} unlabeled, "
                   f"{c['conflict']} conflicting -> {labeled_out}")


@app.command("label-trio")
def label_trio(
    features: Path = typer.Argument(..., help="Feature TSV of the truth-less (wavy) child"),
    parent_gt: list[Path] = typer.Option(
        ..., help="Parent genotype TSVs from `vedette genotype` (repeat for both parents)"
    ),
    out: Path = typer.Option(..., help="Output labeled feature TSV"),
    min_gq: int = typer.Option(20, help="Parental GQ required for carrier/confident-ref"),
    known_sites: Optional[Path] = typer.Option(
        None, "--known-sites",
        help="Reference cohort VCF: calls on artefact-flagged or common sites get no label",
    ),
    max_ref_carrier_frac: float = typer.Option(
        0.5, help="With --known-sites: drop calls whose site is carried by >= this fraction",
    ),
    parent_calls: list[Path] = typer.Option(
        [], help="Each parent's own VCF, in --parent-gt order: TP needs the carrier parent's own call, FP no call in either parent",
    ),
) -> None:
    """Transmission labels from clean parents for a truth-less child (wavy training)."""
    from .train.label_trio import label_by_transmission

    counts = label_by_transmission(
        features, parent_gt, out, min_gq=min_gq,
        known_sites=known_sites, max_ref_carrier_frac=max_ref_carrier_frac,
        parent_calls=parent_calls or None,
    )
    typer.echo(
        f"{counts['tp']} TP, {counts['fp']} FP, {counts['ambiguous']} ambiguous, "
        f"{counts['common']} common/artefact -> {out}"
    )


@app.command()
def train(
    labeled: list[Path] = typer.Argument(..., help="Labeled feature TSV(s) from `vedette label`"),
    out: Path = typer.Option(..., help="Output model directory"),
    model_version: str = typer.Option("v1", help="Version tag written to the manifest"),
    training_data: str = typer.Option("", help="Free-text provenance for the manifest"),
    allow_constant_features: bool = typer.Option(
        False, "--allow-constant-features",
        help="Train even if a feature has no variance (default: refuse — a constant feature usually means broken extraction, as lcr_ovl was through v5)",
    ),
    neutralize: list[str] = typer.Option(
        [], "--neutralize",
        help="NAME=VALUE: pin a feature to a constant at training and scoring (e.g. wave_score=-1; see filter/model.py)",
    ),
) -> None:
    """Train the calibrated DEL/DUP filter models (HistGB + isotonic, chrom-CV)."""
    from .train.train import train_models

    neutral = {}
    for item in neutralize:
        name, _, value = item.partition("=")
        neutral[name.strip()] = float(value)
    manifest = train_models(
        labeled, out, model_version, training_data,
        allow_constant_features=allow_constant_features, neutral_features=neutral,
    )
    for svtype, metrics in manifest["metrics"].items():
        typer.echo(f"{svtype}: PR-AUC(OOF) {metrics['pr_auc_oof']}, n={metrics['n']}")
    typer.echo(f"models -> {out}")


if __name__ == "__main__":
    app()
