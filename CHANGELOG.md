# Changelog

## Unreleased

- **qc** — `vedette wave` writes the wave QC score into evidence stores called
  before schema 0.4 (recomputed from the stored depth, exact to float32). On
  dataset_paper, 2,517 stores had none; 476 of the 2,670 genomes are wavy.
- **report** — from an IGV review of 138 rows: no sub-threshold rows for wavy
  genomes or for older, lower-quality genomes (`--low-quality-prefix`, IP at
  the GHFC: 1/10 and 1/13 true vs 15/28 elsewhere); de novo rows flagged
  « < 1 kb » and « récurrent » below 10 kb (the 7 false, none of the 56 true).
- **report** — the de novo tables list events, not cohort sites: a child's
  sites take its own bounds and nested or bridged sites become one row (a
  PMS proband's 8.2 Mb 22q13 deletion was 61 rows).
- **filter** — production model `models/v5n`: the v5 estimators with the
  per-genome `wave_score` pinned to -1 (`neutral_features`, `vedette train
  --neutralize`). v5 had learned a step on it: genomes with wave >= ~0.12
  passed 10-18k calls instead of ~5k (a third of PMS, 97 dataset_paper
  genomes). Identical to v5 on every other genome.
- **mosaic** — low-fraction clones (~15-25 % of the cells) are no longer
  rejected by a fixed BAF floor: the expected BAF shift is measured on the
  sample's own heterozygous sites, against a depth-aware null (a shallower
  region has a wider BAF spread), and must be consistent window by window.
  On PMS and dataset_paper: no segment lost, 2 new segments in each.
- **report** — the mosaic QC alert gives the segment bounds in Mb with one
  decimal.
- **report** — sub-threshold tier on the family pages: own calls rejected by
  the score alone (LowSQ) with a clear copy number (GQ >= 90, >= 1 kb), hitting
  an established dominant NDD gene and found in no genome outside the family,
  shown greyed as `SOUS-SEUIL NDD` for IGV review (toggle « + sous le seuil »).
  Past 3 such rows a genome is deemed noisy: its rows are masked and its other
  rows carry a QC alert. PMS: 8 rows (3 in SHANK3); dataset_paper: 39 rows in
  32 genomes.

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
