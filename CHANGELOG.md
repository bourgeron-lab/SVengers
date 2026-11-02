# Changelog

## Unreleased

- **report, Nextflow** — a re-filter run (`--steps filter,… --calls_dir`)
  published only the filter outputs under `samples/<S>/`, and the report,
  which reads each evidence store there, silently lost the « Chromosomes
  sexuels » section, the per-genome wave alerts and the depth-based row
  bounds (2670-sample run, 2026-11-02). The pipeline now symlinks each
  sample's evidence store and raw call VCF (+ index) from `--calls_dir` into
  `<outdir>/samples/<S>/` (existing entries are left alone; nothing when
  `--calls_dir` is the run's own `samples/`). `vedette report` warns with the
  count of samples without a store and the sections they miss, and stops
  when no sample has one (`--allow-missing-stores` to report anyway).
- **production** — the Nextflow pipeline calls with a panel of normals: the
  per-bin median of 300 unrelated dataset_paper founders (no SVD). PMS
  re-called with it: cohort sites 104,498 -> 75,335, median apparent de novo
  per child 67 -> 52, Mendelian violations 3.13 % -> 2.61 %, all 79 SHANK3
  families kept, the reviewed 22q13 deletions fully PASS. Common CNVs leave
  the depth channel (they stay detectable through junctions).
- **filter** — large duplications (>= 100 kb, CN >= 3) whose heterozygous
  sites show the allelic imbalance of a copy gain (mean |VAF - 0.5| >= 0.12,
  >= 30 sites, < 50 % segmental duplication) keep PASS whatever their SQ
  (`SQ_OVERRIDE=LARGE_DUP_BAF`). The model scored every large depth-only gain
  near 0: a trisomy 21, two 15q11-q13 duplications, a 7q11.23 duplication, a
  trisomy 13 and a trisomy X were LowSQ. One rescue per genome in median.
- **report** — « Anomalies chromosomiques constitutionnelles »: the mosaic
  scan's whole-copy segments (aneuploidies, arm-level and large gains or
  losses, BAF-confirmed) get a table, and a family-page row when no PASS row
  covers them; fragments fold into it.
- **call** — the wave QC score is measured on the depth as sequenced, before
  any panel of normals (every wave threshold is calibrated on it); with a PoN
  the corrected depth is scored too (`wave_score_pon`). Evidence-store schema
  0.5 (0.1-0.4 still read).
- **report** — row bounds are no longer set on depth for duplications or for
  carriers with wave >= 0.2: a second IGV review of 39 adjusted rows found the
  adjustment right 7/7 on carriers below wave 0.128, 9/18 above 0.2, and both
  adjusted duplications wrong (a coverage wave read as a 23 Mb gain).
- **report** — fragments of a mosaic segment (same lone carrier, inside the
  segment, at its copy number) are folded into the mosaic row; a deeper
  event inside it stays a row. A clonal dataset_paper genome went from 120
  fragment rows to its 12 mosaic rows.
- **cohort stages on several cores** — `merge`, `cohort-filter`,
  `ped-annotate`, `bench trios` and `report` take `--jobs` (default: every
  allowed CPU; Nextflow passes its `cpus`, now 16 for these stages). Each ran
  on one core: ~14 h end to end after the calls on the 2670-sample cohort.
  The work is split by contig, by block of sites, by 10 Mb region or by
  family and restitched in file order, so the outputs are those of the
  serial run (checked on PMS: identical sites, cohort VCF, Mendelian report
  and family pages). PMS with 16 cores: merge 7 min -> 70 s, cohort-filter
  17 min -> 86 s, ped 3.5 min -> 19 s, report 12 min -> 3 min.

- **report** — « Chromosomes sexuels »: per genome, chrX and chrY copy
  number on the raw depth of the bins where the panel of normals tells the
  sexes apart, chrX heterozygosity relative to the autosomes and chrX BAF;
  XX/XY/X0/XXY/XXX/XYY, intermediate levels (mosaic, partial loss of Y), a
  heterozygous single X (DNA contamination?) and pedigree-sex discordances
  are listed and flagged in the carrier QC column. On 2989 GHFC genomes: one
  XXY, one XYY, one XXX, a probable 45,X/46,XX mosaic, 8 partial losses of Y,
  4 men with a heterozygous X, and 5 PMS families whose two parents have
  their sexes swapped.
- **call (fix)** — with a panel of normals, a male chrX/chrY was divided by
  the male panel median (~0.5) and came out at ratio 1.0 against a baseline
  of 0.5: every male had his whole chrX and chrY called a CN 2 DUP (1,501
  dataset_paper men, 156 Mb of chrX in median), and cohort genotyping read
  every man as a carrier of the chrX DUP sites. The normalized ratio is
  brought back to the sample's baseline (ploidy / 2). Re-called male genome:
  autosomes identical, chrX depth DUP 156 -> 0.2 Mb; female genome
  identical. Stores and calls made with a PoN before this fix must be
  re-called for the men.
- **call, genotype (fix)** — the depth noise (sigma0) is fitted around the
  contig's own baseline and expressed at CN2: on a male chrX at 0.5, the
  former 0.7-1.3 window kept ~170 of 147k bins and gave sigma0 0.02 instead
  of 0.08, so single-bin dips became hemizygous deletions (about a third of
  the men genotyped 1/1 at recurrent chrX sites, ~75 chrX "de novo" per
  child). Hidden before the PoN fix, which left male chrX at 1.0.
- **report** — the de novo candidates a sex-chromosome aneuploidy explains
  (chrX of an XXY, XXX or X0 child, chrY of an XYY one) leave the de novo
  guard's count and the family tables, with a note pointing to « Chromosomes
  sexuels ». An XXY boy of dataset_paper had 681 of his 695 de novo sites on
  chrX: he was excluded as a technical outlier, and a 1.4 Mb piece of his
  extra X over MECP2 read as a de novo duplication. On the family events
  pages those calls read « caryotype XXY » instead of « de novo », without
  the de novo alerts.
- **karyotype** — GC bias guard: raw depth carries the genome's GC bias and
  chrX, AT-rich, reads like chr4. A dataset_paper man (wave 0.50, autosomes
  from 1.68 copies on chr4 to 3.55 on chr19) read chrX 0.84 / chrY 0.65
  (« intermédiaire »), a PMS woman chrX 1.70 (« mosaïque ? »). New
  `gc_contrast` (median CN of chr16/17/19/22 minus chr4/5/13/18; 14 of 2989
  genomes above 0.5): above 0.5 an abnormal X/Y reading becomes « ? » with a
  GC-bias note, an XX/XY one is kept. Only those two genomes change.
- **cohort-filter, ped-annotate, bench trios** — chrX/chrY aware: men are
  hemizygous outside the PARs (male founders out of the HWE test, carrier
  fraction per sex, CN discordance against the carrier's ploidy), a son's
  chrX comes from his mother only, chrY from the father only.
- **cohort-filter** — `ChrYNoDepth`: on chrY outside the PARs, a site whose
  carriers' copy number does not support the SV type (fewer than half of
  them above / below their baseline) leaves PASS. dataset_paper, once the
  men's chrY was called at its baseline: 525 chrY PASS sites, 453 of them
  small junction-only calls in the Y repeats carried by up to two thirds of
  the men, with a son disagreeing with his father in 40-54 % of the pairs
  carrying them; the 70 deletions the copy number supports, 20 %. Sites
  shorter than one depth bin (1 kb) are flagged too, their carriers' CN
  being set by the junctions: father/son discordance 34 % (DEL) and 86 %
  (DUP) below 1 kb, 10 % and 24 % above; 55 chrY PASS sites left.
- **filter** — the default operating threshold is SQ >= 5 (was 13), the
  operating point of the shipped model `models/v5n`: without a site TOML
  setting it, a run filtered far more than production.
- **pedigree** — sex written as "1.0"/"2.0" (the GHFC pedigrees) is parsed;
  it read as unknown for every individual.

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
