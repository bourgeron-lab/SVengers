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
def merge(
    vcfs: list[Path] = typer.Argument(..., help="Filtered single-sample VCFs (or one @filelist)"),
    out: Path = typer.Option(..., help="Output sites-only VCF (.vcf.gz)"),
    fasta: Path = typer.Option(..., help="Reference FASTA (for the contig header)"),
    min_gq: int = typer.Option(0, help="Minimum GQ for a call to nominate a site"),
    include_filtered: bool = typer.Option(
        False, help="Also cluster non-PASS calls (default: PASS/LowPosterior only)"
    ),
) -> None:
    """Cross-sample clustering of filtered calls into cohort sites."""
    from .cohort.merge import run_merge
    from .reference import read_fai

    paths: list[Path] = []
    for v in vcfs:
        if str(v).startswith("@"):
            paths.extend(Path(p) for p in Path(str(v)[1:]).read_text().split())
        else:
            paths.append(v)
    n = run_merge(
        paths, read_fai(fasta), out, min_gq=min_gq, pass_only=not include_filtered
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
) -> None:
    """Assemble genotype tables into a cohort VCF with AF/HWE/artifact metrics."""
    from .cohort.metrics import run_cohort_filter
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
    )
    typer.echo(
        f"{counts['sites']} sites, {counts['artifact']} CohortArtifact -> {out}"
    )


@app.command("ped-annotate")
def ped_annotate(
    cohort_vcf: Path = typer.Argument(..., help="Multi-sample cohort VCF"),
    pedigree: Path = typer.Option(..., help="Pedigree TSV"),
    out_report: Path = typer.Option(..., help="Per-trio-site report TSV"),
    out_summary: Optional[Path] = typer.Option(None, help="Per-trio summary TSV"),
) -> None:
    """Mendelian consistency, inheritance and de novo annotation on a cohort VCF."""
    from .cohort.mendel import run_ped_annotate
    from .ped import parse_ped

    stats = run_ped_annotate(
        cohort_vcf, parse_ped(pedigree), out_report, out_summary
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


@bench_app.command("trios")
def bench_trios(
    cohort_vcf: Path = typer.Argument(..., help="Multi-sample cohort VCF"),
    pedigree: Path = typer.Option(...),
    out: Path = typer.Option(..., help="Output sweep TSV"),
    field: str = typer.Option("GQ", help="FORMAT field to threshold (GQ or SQ)"),
    plot: bool = typer.Option(False, help="Also write MVR/transmission figure (PNG)"),
) -> None:
    """MVR and transmission rate as a function of the quality threshold."""
    from .bench.trios import sweep
    from .ped import parse_ped

    rows = sweep(cohort_vcf, parse_ped(pedigree), field=field, out_tsv=out)
    if plot:
        from .bench.plots import plot_trio_sweep

        plot_trio_sweep(rows, Path(out).with_suffix(".png"))
    typer.echo(f"sweep -> {out}")


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
def train(
    labeled: list[Path] = typer.Argument(..., help="Labeled feature TSV(s) from `vedette label`"),
    out: Path = typer.Option(..., help="Output model directory"),
    model_version: str = typer.Option("v1", help="Version tag written to the manifest"),
    training_data: str = typer.Option("", help="Free-text provenance for the manifest"),
) -> None:
    """Train the calibrated DEL/DUP filter models (HistGB + isotonic, chrom-CV)."""
    from .train.train import train_models

    manifest = train_models(labeled, out, model_version, training_data)
    for svtype, metrics in manifest["metrics"].items():
        typer.echo(f"{svtype}: PR-AUC(OOF) {metrics['pr_auc_oof']}, n={metrics['n']}")
    typer.echo(f"models -> {out}")


if __name__ == "__main__":
    app()
