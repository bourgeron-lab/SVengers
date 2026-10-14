# Changelog

## Unreleased

- **report** — row bounds are no longer set on depth for duplications or for
  carriers with wave >= 0.2: a second IGV review of 39 adjusted rows found the
  adjustment right 7/7 on carriers below wave 0.128, 9/18 above 0.2, and both
  adjusted duplications wrong (a coverage wave read as a 23 Mb gain).
- **report** — fragments of a mosaic segment (same lone carrier, inside the
  segment, at its copy number) are folded into the mosaic row; a deeper
  event inside it stays a row. A clonal dataset_paper genome went from 120
  fragment rows to its 12 mosaic rows.

## v1.1.0 — 2026-10-13

Production model, mosaic detection and a report reworked after an IGV review
of 138 rows by the GHFC lab. Same evidence-store and feature schemas as
v1.0.0: existing stores and calls stay valid.

### Filtering

- **Production model `models/v5n`**: the v5 estimators with the per-genome
  `wave_score` pinned to -1. v5 had learned a step on that feature: genomes
  with wave >= ~0.12 passed 10-18k calls instead of ~5k (a third of PMS, 97
  dataset_paper genomes). Identical to v5 on every other genome. On PMS:
  cohort sites 215,867 -> 104,498, median apparent de novo per child
  273 -> 67; on dataset_paper: 522,992 -> 361,359 sites, 559 -> 194.
- `neutral_features` in the model manifest and `vedette train --neutralize
  NAME=VALUE` pin a feature to a constant at training and scoring.

### Mosaics

- Low-fraction clones (~15-25 % of the cells) are no longer rejected by a
  fixed BAF floor: the expected BAF shift is measured on the sample's own
  heterozygous sites, against a depth-aware null (a shallower region has a
  wider BAF spread), and must hold window by window. PMS and dataset_paper:
  no segment lost, 2 new in each (a 20 % 22q13 loss confirmed in IGV).

### QC

- `vedette wave` writes the coverage-wave score into evidence stores called
  before schema 0.4, recomputed from the stored depth (exact to float32).
  dataset_paper: 2,517 stores had none; 476 of its 2,670 genomes are wavy.

### Report

- **Row bounds set on the carriers' depth** for events of 50 kb and more,
  compared with the family's non-carriers: on the reviewed rows, the start
  error went from 0.3-10 Mb to a few kb (e.g. 50,435,001 -> 49,884,001 for
  49,883,222 in IGV).
- **De novo tables list events, not cohort sites**: a PMS proband's 8.2 Mb
  22q13 deletion was 61 rows.
- **De novo alerts** « < 1 kb » and « récurrent » below 10 kb: the 7 false
  high-impact de novo rows of the review, none of the 56 true ones.
- **Sub-threshold NDD tier**: calls rejected by the score alone with a clear
  copy number in an established NDD gene, private to the family, shown
  greyed for review; none for wavy genomes or for older, lower-quality
  genomes (`--low-quality-prefix`), at most 3 per genome.
- Mosaic segments are rows of the family pages; the mosaic alert gives the
  bounds in Mb with one decimal.
- No dead link for families above the per-page member limit.
- Lower memory: sub-threshold candidates held in compact arrays (dataset_paper
  went past 32 GB).

## v1.0.0 — 2026-10-05

First release. Germline DEL/DUP calling from WGS with post-calling filtering,
cohort re-genotyping and pedigree-aware reporting, run end to end by the
Nextflow wrapper (`workflows/nextflow`) on SLURM.

### Pipeline

- **call** — per-sample multi-evidence calling (depth CN-HMM, discordant
  pairs, split reads incl. hidden splits, insert-size shift, BAF) into a
  schema-versioned Zarr evidence store; downstream stages never re-read the
  CRAM.
- **filter** — hard filters, then a calibrated per-type score (SQ,
  gradient boosting, models/v5, operating threshold SQ ≥ 5); large deletions
  the model under-scores are kept when chromosome-level depth agrees.
- **merge / genotype** — cohort site clustering with recurrence rescue, then
  re-genotyping of every sample from its store.
- **cohort-filter** — artefact, HWE, carrier-fraction and CN-discordance
  flags; directional genotype discordance; reference-cohort annotation
  (`--known-sites`, CommonRef).
- **ped** — Mendelian inheritance and de novo candidates per trio.
- **mosaic** — chromosomal-scale mosaic CNVs (fractional copy number,
  confirmed by BAF and window homogeneity), per sample from the store.
- **report** — cohort HTML report and per-family event pages.

### Family event pages

- One row per biological event: overlapping/adjacent sites of the same
  carriers merged, ends taken from the carriers' own chained PASS calls,
  genes/impact recomputed, majority inheritance.
- Impact score after the ACMG/ClinGen CNV grid (sections 1-3), DEL and DUP
  scored separately, GeneTrek standing in for the ClinGen dosage map.
- Carrier QC alerts: wave (without the row's own contig), PASS carrier load
  (≥ 5× cohort median), de novo outlier (≥ 5× cohort median), large-CNV
  burden (≥ 3 events ≥ 1 Mb), mosaic segment.
- Own-call status, re-genotyped CN support, LCR/segdup overlap, IGV links,
  manual true/false/uncertain validation with TSV export/import.

### Known limitations

- Specificity of carriers assigned by cohort re-genotyping alone (no own
  call) is low; the pages show the own-call status and filter on it by
  default.
- A carrier without an own PASS call keeps its cohort site's breakpoints,
  which can overshoot its event (edge-CN test at genotyping planned).
- Mosaic detection: autosomes only, segments ≥ 3 Mb; complete losses are not
  confirmed by BAF (no heterozygous sites left).
- Impact score: no triplosensitivity (gains of dosage-sensitive genes are
  under-scored), no 5'/3' distinction for partial deletions, inheritance and
  phenotype not counted.
- The stored wave score is inflated by large CNVs; the report corrects it
  per row only.
- Per-sample filtered VCFs written before v1 may carry call ends up to 1 kb
  past a contig end; merge and report bound them.
- Not published to PyPI or bioconda.
