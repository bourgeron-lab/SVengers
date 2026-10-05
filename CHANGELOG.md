# Changelog

## v1.1.0 — 2026-10

### Filter model v6 (production)

- `models/v6` replaces `models/v5` (operating threshold unchanged, SQ ≥ 5).
  Same training data and labels; the GIAB samples are re-called with the
  current code and features recomputed. v5 had been trained on GIAB calls
  made before the negative-binomial junction likelihood changed `gq`, so in
  production it served a feature whose distribution had moved.
- Holdout of 8 wavy GHFC children (production features, unseen by both
  models), SQ ≥ 5: DEL precision/recall 0.974/0.714 → 0.978/0.973, DUP
  0.795/0.728 → 0.878/0.960.
- HG002, first out-of-sample figure (v6 trained without HG002): DEL
  precision 0.870, recall 0.365, F1 0.514 (GIAB v5.0q); best F1 of the
  compared callers from 300 bp up; CMRG F1 0.573. Earlier HG002 figures were
  in-sample. Tables: `docs/benchmark_hg002_v6.tsv`.

### Fixes

- `bench`: truth VCFs that declare no SVTYPE/SVLEN (CMRG) are read.
- `configs/giab.toml` carries the production operating point; GIAB runs keep
  1 kb depth bins (a CRAM-only run with `ghfc.toml` scans at 100 bp).

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
  gradient boosting, models/v5 in v1.0.0, operating threshold SQ ≥ 5); large deletions
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
