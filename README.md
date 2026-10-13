# Vedette

**Vedette** is a germline copy-number variant caller (deletions and duplications) for
short-read whole-genome sequencing, built for **specificity** at cohort scale.

It combines read depth, discordant pairs, split reads, insert-size shifts and B-allele
frequency per sample, filters the calls with a calibrated machine-learning score, then
re-genotypes every cohort site in every sample and reads the result through the pedigree:
inheritance, de novo candidates, mosaic chromosomal segments, and an HTML report with one
interactive page per family.

- **One CRAM pass per sample.** Calling writes a compact, schema-versioned
  [Zarr](https://zarr.dev) *evidence store* (~100 MB). Filtering, cohort re-genotyping,
  mosaic scanning and QC all read the store and never touch the CRAM again.
- **Tiered filtering.** Hard filters, then a calibrated gradient-boosting score (`SQ`,
  phred-scaled probability that the call is real), then cohort metrics (recurrence,
  Hardy-Weinberg, carrier fraction, copy-number consistency).
- **Built for families.** Mendelian consistency, transmission and de novo quality on
  trios, with guards against technically atypical genomes (carrier load, de novo excess,
  wavy coverage).
- **Scales from one genome to thousands**, through a Nextflow wrapper on SLURM.

Status: **v1.1.0**, in production at the GHFC lab (Institut Pasteur) on cohorts of
300 to 2,700 genomes. Not yet on PyPI or Bioconda; install from source (below).

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [How it works](#how-it-works)
- [Commands](#commands)
- [Filter models and benchmark](#filter-models-and-benchmark)
- [Outputs and the HTML report](#outputs-and-the-html-report)
- [Configuration](#configuration)
- [Development](#development)
- [Citation and license](#citation-and-license)

## Installation

Vedette needs Python ≥ 3.11 and a few scientific packages (pysam, numpy, scipy, numba,
zarr, scikit-learn). The conda/micromamba environment in `containers/environment.yml`
pins them:

```bash
git clone https://github.com/darkxander/vedette.git
cd vedette
micromamba create -n vedette -f containers/environment.yml
micromamba run -n vedette pip install -e .
micromamba run -n vedette vedette --help
```

A container is also provided: `containers/Dockerfile` (Docker) and `containers/vedette.def`
(Apptainer/Singularity).

Training and benchmarking (not needed to call) additionally use
[Truvari](https://github.com/ACEnglish/truvari): `pip install truvari`.

## Quick start

### One sample

```bash
# 1. call: CRAM -> raw calls (VCF) + evidence store
vedette call --sample S1 --cram S1.cram --fasta GRCh38.fa \
    --config my_site.toml --resources resources/ --out-dir out/S1

# 2. filter: hard filters + calibrated SQ score (PASS = SQ >= 5 with the shipped model)
vedette filter out/S1/S1.vedette.vcf.gz --store out/S1/S1.evidence.zarr \
    --model-dir models/v5n --config my_site.toml --out out/S1/S1.filtered.vcf.gz

# 3. optional: chromosomal-scale mosaic segments, from the store alone
vedette mosaic out/S1/S1.evidence.zarr --out out/S1/S1.mosaic.tsv
```

If your pipeline already produces a 1 kb depth track and a per-site VAF track, pass them
with `--depth-bedgraph` / `--vaf-bedgraph` alongside the CRAM: depth and BAF are then read
from the tracks (the "fast path") and the CRAM is scanned for junction evidence only. Our
production tracks are `mosdepth -n -m -x -Q 40 --by 1000` (median depth per 1 kb bin) and
DeepVariant 1.9 sites with `FILTER=PASS`, `QUAL ≥ 30`, `DP ≥ 10`.

### A cohort

The Nextflow wrapper in `workflows/nextflow/` runs every step on SLURM, from a sample
manifest (`vedette manifest discover` builds one) and a 6-column pedigree:

```bash
cd workflows/nextflow
./launch_driver.sh --manifest cohort.manifest.tsv --pedigree cohort.pedigree.tsv \
    --outdir results/cohort --cohort_name "My cohort" \
    --steps call,filter,merge,genotype,cohort,ped,report \
    --family_events true --events_pass_only true
```

Steps can restart at `filter` from published calls (`--calls_dir`), e.g. after a model
update. Set the site-specific paths (reference, masks, model, resources) in
`nextflow.config`.

## How it works

```
          per sample                                 cohort
 CRAM ──► call ──► evidence store ──► filter ──► merge ──► genotype ──► cohort-filter ──► ped-annotate ──► report
 (+ optional        (Zarr: depth, BAF,   (hard filters   (cluster    (every site in     (AF, HWE,         (Mendel,         (HTML,
  depth/VAF          pairs, splits,       + SQ score)     PASS calls  every sample,      artefacts,        de novo,         family pages)
  tracks)            insert sizes)                        into sites) from its store)   carrier frac.)    inheritance)
                          │
                          └──► mosaic (fractional copy number confirmed by BAF)
```

- **Calling.** Depth is GC-corrected and normalised to the diploid autosomal median (male
  chrX at 0.5); an integer copy-number HMM segments it. Junction evidence (discordant pairs,
  split reads, including hidden splits in soft clips) nominates breakpoints. Calls stay
  permissive at this stage; specificity is the filter's job.
- **Filtering.** Hard filters (masks, mappability, missing evidence), then the calibrated
  `SQ` score from 30 features. The feature schema is append-only and each model pins its
  hash.
- **Cohort.** Filtered calls are clustered into sites; every site is re-genotyped in every
  sample from its evidence store, so carriers missed at the single-sample stage are
  recovered without re-reading any CRAM.
- **Pedigree.** Mendelian consistency, transmission, de novo quality (`DNQ`), plus guards
  for technically atypical genomes.
- **Mosaics.** A 1 Mb-window scan finds chromosomal-scale segments at fractional copy
  number (clones present in 15 to 85 % of cells) and keeps them only when the allelic
  imbalance matches the fraction, against a depth-aware null (GC waves widen the BAF
  spread without any imbalance).

## Commands

| command | role |
|---|---|
| `vedette call` | per-sample calling: CRAM (and optional tracks) → raw VCF + evidence store |
| `vedette filter` | hard filters + calibrated `SQ` score → filtered VCF (and feature matrix) |
| `vedette mosaic` | chromosomal-scale mosaic segments from a store, or a whole run (`--run-dir`) |
| `vedette wave` | write the coverage-wave QC score into stores that lack it (pre-0.4 stores) |
| `vedette merge` | cluster filtered calls across samples into cohort sites |
| `vedette genotype` | re-genotype every cohort site in one sample, from its evidence store |
| `vedette cohort-filter` | cohort VCF with allele frequency, HWE, artefact and consistency metrics |
| `vedette ped-annotate` | Mendelian consistency, inheritance, de novo annotation |
| `vedette ped` | parse and summarise a pedigree |
| `vedette report` | self-contained HTML cohort report and per-family pages |
| `vedette manifest` | build (`discover`) or check (`validate`) a cohort sample manifest |
| `vedette pon` | panel of normals: depth profile and recurrent-artefact catalogue |
| `vedette label` / `label-trio` | training labels: Truvari against a truth set, or transmission in trios |
| `vedette train` | train the calibrated DEL/DUP models (`--neutralize` pins a feature) |
| `vedette bench` | benchmark harness: truth-set evaluation, caller comparison, trios |

Every command has `--help`.

## Filter models and benchmark

The production model is **`models/v5n`**: gradient boosting with isotonic calibration,
trained on GIAB HG002/HG003/HG004 (labels from the GIAB v5.0q truth with Truvari, and by
transmission for the parents) and on 16 GHFC children with wavy coverage (transmission
labels). PASS is `SQ ≥ 5`.

`v5n` is v5 with the per-genome `wave_score` feature neutralised: v5 had learned a step on
it and passed 10 to 18 thousand calls on genomes with strong coverage waves instead of
about 5 thousand. On every other genome, the two are identical. See `docs/ROADMAP.md` for
the history of models v1 to v8 and why the later ones were not retained.

**HG002 (GIAB v5.0q stvar and CMRG v1.00), deletions, Truvari:**

| truth | caller | recall | precision | F1 |
|---|---|---|---|---|
| GIAB v5.0q | **Vedette (v5, SQ ≥ 5)** | 0.334 | **0.951** | 0.495 |
| GIAB v5.0q | Manta | 0.376 | 0.934 | 0.536 |
| GIAB v5.0q | dysgu | 0.389 | 0.879 | 0.539 |
| GIAB v5.0q | DELLY | 0.284 | 0.826 | 0.423 |
| CMRG v1.00 | **Vedette (v5, SQ ≥ 5)** | 0.370 | 0.944 | 0.531 |
| CMRG v1.00 | Manta | 0.489 | 0.978 | 0.652 |
| CMRG v1.00 | dysgu | 0.457 | 0.840 | 0.592 |
| CMRG v1.00 | DELLY | 0.261 | 0.857 | 0.400 |

HG002 is part of the training set, so these figures are in-sample; out-of-sample
estimates for retrained models are in `docs/ROADMAP.md`. Per-size and per-context
tables: `docs/benchmark_hg002_v1.tsv`. Vedette's strength is precision from 300 bp to
100 kb (F1 0.72 at 300 bp–1 kb, 0.84 at 1–10 kb). Its weak spots are small deletions
(50–300 bp) and segmental duplications. On real cohorts the decisive metrics are trio-based
(apparent de novo rate, Mendelian violations), which the report shows for every family.

## Outputs and the HTML report

A cohort run writes, under `--outdir`:

```
samples/<S>/<S>.vedette.vcf.gz      raw calls
samples/<S>/<S>.evidence.zarr       evidence store
samples/<S>/<S>.filtered.vcf.gz     filtered calls (FILTER, SQ)
samples/<S>/<S>.mosaic.tsv          mosaic segments
cohort/cohort.vcf.gz                re-genotyped cohort VCF with metrics
cohort/mendel_report.tsv            per-trio, per-site inheritance and de novo
report.html                         cohort report
families/<FID>/events.html          one interactive page per family
```

The **family pages** list every event carried by a member, with:
- bounds taken from the carrier's own calls and, from 50 kb, set on the carriers'
  depth against the rest of the family, so they match IGV;
- genes and an ACMG/ClinGen-style **impact score** (dosage-sensitive and established
  neurodevelopmental genes, gene count);
- inheritance, whether the member has an own call or only a re-genotyped one, and carrier
  QC alerts (wave, carrier load, de novo excess, large-CNV burden, mosaic segment,
  small or recurrent de novo);
- rows for mosaic segments, and an optional review tier of sub-threshold calls in
  established NDD genes;
- filters, an IGV button (IGV desktop batch port) and vrai/faux/incertain validation boxes
  exportable as TSV.

## Configuration

Site-specific resources live in a TOML file passed with `--config` (reference FASTA,
exclusion, segmental-duplication and low-complexity masks, population-coverage track,
filter threshold). `configs/ghfc.toml` and `configs/giab.toml` are complete examples to
adapt; the Nextflow parameters are in `workflows/nextflow/nextflow.config`.

## Development

```bash
micromamba run -n vedette python -m pytest tests/ -q
micromamba run -n vedette ruff check src tests
```

Design notes, the decisions behind each component and the open work are in
`docs/ROADMAP.md`; release notes are in `CHANGELOG.md`.

Invariants worth knowing before contributing:
- downstream stages never re-read a CRAM;
- any change to the store's arrays or attributes bumps its schema version;
- the feature schema is append-only;
- positions are 0-based half-open internally.

## Citation and license

If you use Vedette, please cite it with the metadata in `CITATION.cff`.

MIT License — Alexandre Mathieu, Institut Pasteur (GHFC lab).
