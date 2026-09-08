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
        model_dir=model_dir, features_out=features_out,
    )
    typer.echo(
        f"{counts['total']} calls: {counts['pass']} PASS, "
        f"{counts['hard_filtered']} hard-filtered, {counts['low_sq']} LowSQ -> {out}"
    )


@app.command()
def merge() -> None:
    """Cross-sample clustering of filtered calls into cohort sites."""
    _not_yet("merge")


@app.command()
def genotype() -> None:
    """Force-genotype cohort sites in every sample from evidence stores."""
    _not_yet("genotype")


@app.command("cohort-filter")
def cohort_filter() -> None:
    """Cohort-level AF/HWE/artifact/batch annotation."""
    _not_yet("cohort-filter")


@app.command()
def pon() -> None:
    """Panel-of-normals depth profile and artifact catalogue."""
    _not_yet("pon")


if __name__ == "__main__":
    app()
