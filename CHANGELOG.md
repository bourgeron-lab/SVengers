# Changelog

## Unreleased

- **`cohort-filter`: a small cohort takes the reference cohort's
  CohortArtifact** — with `--known-sites`, a site common in the cohort (AF >
  `--artifact-af`) whose matching reference site is CohortArtifact is
  CohortArtifact too (INFO `REF_ARTIFACT`; `--no-known-artifacts` turns it
  off). Rare sites keep their own verdict. The founders' HWE test lacks power
  and stability on a small cohort: leaving out any one founder flips ~150
  PMS sites (206 founders), 59 on dataset_paper (1,688), and of PMS's 54,492
  common sites matched in dataset_paper, dataset_paper flags 10,767 that PMS
  does not. On PMS: 11,597 sites flagged from the reference, PASS sites
  28,093 -> 23,129 (all lost ones common: carrier fraction 0.03-0.6, median
  0.33), family-page rows -39 %, PASS de novo candidates 2,203 -> 1,302; no
  22q13 deletion touched (lost chr22 sites > 45 Mb all < 2 kb).

- **`vedette regenotype`** — recomputes GT, GQ, CN, CNQ, RD, BAFN and BAFDEV
  of a sample's raw calls from its evidence store, junction support and
  source kept (INSERT calls copied), so that the filter can be re-run on a
  BAF rewritten by `vedette baf`. On stores whose BAF was not rewritten it
  reproduces the calling-time values: GT, CN and BAFN identical on 44,895
  calls of two genomes, GQ on all but 2, and the filtered output is the same.

## v1.9.0 — 2026-11-05

- **`vedette labels --from-genetics-viz`** (#11) — reads the curators'
  validations of genetics-viz (`validations/vedette.tsv`, append-only: the
  last row of a family, variant and sample not marked Ignore wins) into one
  label table (`present` = 1, `absent` = 0; `uncertain` and `different`
  listed, no training label), with the curated-minus-called breakpoint
  offsets; `--model-version` keeps the validations of one model, and
  `--features` labels one sample's feature rows for `vedette train` (same
  type, reciprocal overlap >= 0.5).
- **Family export: per-carrier evidence** (#8) — for each member, the own
  call behind the row (its filtered VCF, same type, reciprocal overlap >=
  0.5, a PASS call first): `OWN_`, `SQ_`, `SRC_`, `PRECISE_` columns in the
  TSV, and in the Parquet file its id, FILTER, bounds, CIPOS/CIEND, GQ, PE,
  SR, CG, RD, CNQ and model version; per row, the cohort site's AF,
  NSAMPLES, PRECISE_FRAC, SQMEAN, HWE_P, CN_DISC and REF_AF.
- **Site keys** (#9) — `site_key` (`chr:start-end:del|dup`, the cohort
  site's own 1-based bounds) in the family export, and `vedette site-map OLD
  NEW` for the correspondence of sites across a re-merge: same bounds first,
  then mutual best matches by reciprocal overlap (>= 0.5), then leftovers as
  split or merged; on two identical cohort VCFs every site is `same`
  (pairing every overlap read 92,004 PMS sites as "complex": sites overlap
  their neighbours in a real cohort).
- **Family export as Parquet** (#6) — `<FID>.vedette.parquet` with pyarrow
  (`vedette[parquet]`): one line per row x member, typed, list columns,
  versions and the panel of normals in the file metadata; the `.json` gives
  `pon` too.

## v1.8.0 — 2026-11-05

- **`vedette baf`: BAF without a QUAL threshold** — rewrites the BAF of
  existing evidence stores from the sample's SNV VCF (DeepVariant: PASS
  biallelic SNVs, DP >= 10, no QUAL threshold), without the CRAM, about a
  minute per genome; the previous arrays stay as `<contig>/baf_track/`. The
  GHFC VAF track keeps QUAL >= 30, which drops many true heterozygotes: on
  copies of production stores the autosomal het/hom ratio goes from 0.04 to
  0.91 for an old (IP) genome, 0.08 to 1.18 for a 91x one, 1.40 to 1.70 for
  ordinary genomes; karyotypes are unchanged, and a man's X heterozygosity
  falls from 0.08-0.12 to 0.02-0.03 on the old genomes. Store schema 0.6
  (`baf_source`). The calls' BAF features (BAFN, BAFDEV) are unchanged
  until a new model is trained on them.
- `--snv-vcf` reads PASS (or unfiltered) sites only, with FORMAT/DP as the
  site depth when given: DeepVariant's RefCall records are no heterozygous
  sites.

## v1.7.0 — 2026-11-05

- **Family export: rare non-PASS sites** (#10) — with `--events-pass-only`
  the events page keeps PASS sites only, and `<FID>.vedette.tsv` now also
  holds the non-PASS sites carried by at most 20 % of the cohort
  (`--export-nonpass-max-carrier-frac`, recorded in the `.json`), for a
  reader to show behind a toggle. Most non-PASS sites are carried by over
  half the cohort (dataset_paper: 167,560 of 190,367, 33x the carriers of the
  PASS sites) and stay out; the export grows by about 7 %.
- **Report: one row per inherited sub-threshold event** — members of a
  family carrying the same sub-threshold event (same type, overlap > 50 % of
  the larger) share one row instead of one row each (in v1.6.0, 7 PMS
  events took 15 of the 41 large-DUP rows).

## v1.6.0 — 2026-11-05

- **Report: large sub-threshold duplications** — the family pages keep, as
  "SOUS-SEUIL GRANDE DUP" rows, a carrier's LowSQ duplications (CN >= 3,
  GQ >= 90) once joined across gaps of at most 10 % of the joined span
  (1 Mb at most), when >= 100 kb, private to the family (the joined event
  and its largest piece), on a genome neither wavy nor low quality, and not
  on a sex chromosome the carrier's karyotype changes. The v5 IGV review
  found 12 of 15 sub-threshold DUPs >= 100 kb true, and SQ ranking nothing
  below the threshold; a de novo 15q11.2q13 duplication at CN 4 was cut
  into LowSQ pieces (1.1 and 5.7 Mb). Same per-genome cap as the NDD tier.

## v1.5.0 — 2026-11-05

- **Impact score: ClinGen dosage map** — `vedette report --clingen-regions
  --clingen-genes` read the ClinGen curation lists (GRCh38). A CNV covering
  >= 90 % of a region scored 3 (HI for a loss, TS for a gain) gets +1.00
  (2A), >= 20 % +0.45 (partial overlap kept visible: depth calls stop in the
  repeats flanking recurrent regions); one lying inside a region curated
  "dosage sensitivity unlikely" (score 40) gets -1.00, new level `bénin` /
  `benign`. Gene HI 3 counts as an established dominant NDD gene for losses,
  gene TS 3 gives +1.00 to a complete duplication. GeneTrek still scores the
  genes ClinGen has not. Found by the v5 IGV review: a de novo 15q11.2q13
  duplication (copy number 4, affected child) was on its family page at
  +0.15 ("TS non évaluée"); it is now +0.45, a complete 16p11.2 BP4-BP5
  duplication +1.00 instead of +0.15. On PMS's PASS sites: 14 levels up,
  212 to benign, none down otherwise. The export JSON names the ClinGen files.

## v1.4.2 — 2026-11-05

- **docs** — ROADMAP: `CohortArtifact` is unstable on a small cohort. On
  PMS (214 founders), about 140 common sites sit right at the HWE threshold
  (1e-6), so one founder fewer flips 163 FILTERs; rare sites are unaffected.
  Options listed: a one-sided heterozygote-excess test, a threshold scaled to
  the founder count, an F margin. No change to the package.

## v1.4.1 — 2026-11-05

- **docs** — `docs/EXPORT.md` gave the `wavy` QC flag at wave ≥ 0.2; the
  report flags it at ≥ 0.128 (`qc.WAVE_WARN`). No change to the package.

## v1.4.0 — 2026-11-04

VCF coordinates follow the spec, and the family pages get a machine-readable
export (genetics-viz issues #5, #7, #10). Same store schema, feature schema
and production model (v5n) as v1.3.3.

- **VCF coordinates (format change)** — `POS` is now the padding base before
  the event and `END` its last base (VCF 4.3, symbolic alleles), with
  `SVLEN = ±(END - POS)`. Vedette wrote `POS` on the first affected base and,
  htslib setting `END = POS + |SVLEN|`, `END` one base past the last one: one
  base to the right on both sides (a junction deletion's `POS` sat one base
  after the breakpoint read from the clipped reads). The merge also read the
  per-sample `END` as the end and wrote it +1 again, so a cohort site ran one
  base past its calls. All reading and writing of SV coordinates now goes
  through `vcfout.span` / `vcfout.sv_record`. New files carry
  `##vedette_coords=padded`; Vedette files without it are read with the old
  convention, so existing per-sample VCFs and cohorts stay usable. TSV outputs
  keep their 0-based half-open intervals; their `end` loses the extra base.
- **Family export** — `vedette report --family-events` writes
  `families/<FID>/<FID>.vedette.tsv` (the events page's rows: `kind` site /
  subthreshold / mosaic / chromosomal, cohort site ids, 1-based inclusive
  bounds and `vcf_pos`, `impact_score` / `impact_level` / `impact_reason`,
  complete exonic and intronic coding-gene lists with Ensembl ids, own-call
  and support codes, per-member `GT`/`GQ`/`CN`, each child's inheritance code
  and DNQ), `<FID>.vedette.qc.tsv` (per member: wave, PASS sites carried,
  carrier excess, de novo candidates, karyotype, flag codes) and
  `<FID>.vedette.json` (`schema_version` 1, Vedette and model versions,
  member order). Format and coordinate convention: `docs/EXPORT.md`.
- **Report** — merged rows keep every cohort site id through the second
  collapse (the count of collapsed sites was reset to the rows of the second
  pass).

## v1.3.3 — 2026-11-04

- **pedigree** — a repeated sample line whose sex or phenotype differs from
  the first is no longer reported as an "identical" duplicate: the warning
  names the conflict and the line kept. In the PMS pedigree, two
  individuals were each written twice, once per sex — the second line was a
  missing parent's. The first line is still the one used.

## v1.3.2 — 2026-11-04

- **tests** — `tests/test_tensor.py` passes the CI lint (import order,
  `zip(strict=True)`); v1.3.0 and v1.3.1 failed the lint job, tests and
  container passed. No change to the package.

## v1.3.1 — 2026-11-04

- **tensors (fix)** — counts are capped at 60,000 before the float16 cast:
  soft-clip piles in repeats passed float16's 65,504 and were written as
  inf (15 of the first 299 production genomes, `sr_left`).

## v1.3.0 — 2026-11-04

First step of the CNN second opinion (ROADMAP 1sexdecies). No change to
calling, filtering or the cohort stages; same store, feature schema and
production model as v1.2.0.

- **tensors** — new `vedette tensors`: each call of a filtered VCF (default
  PASS and LowSQ, >= 1 kb) becomes a 13 x 256 track over the call and one
  call length of flank per side (>= 5 kb), read from the evidence store
  only — log2 copy ratio and masked fraction, the panel of normals' median
  and MAD on the same bins (cohort context), heterozygous-site density and
  BAF deviation, DEL/DUP pair ends, soft clips by side, CIGAR deletions,
  insert-size shift, call footprint — plus 12 scalars. One `.npz` per
  genome (~1 MB, ~20 s, ~1,200 calls). The PMS SHANK3 deletion rejected by
  v5n (SQ 0.38) reads as a clean two-step loss: log2 ratio -1.0 inside,
  flat flanks, no heterozygous site, flat panel.

## v1.2.0 — 2026-11-04

Panel of normals in production, sex chromosomes handled end to end, and a
sex-chromosome karyotype per genome. dataset_paper (2670 genomes) and PMS
(319) re-called and re-run with this version: autosomes unchanged against
the previous production, chrX of men at their hemizygous baseline,
Mendelian violation rate 3.15 % -> 2.38 % on dataset_paper. Evidence-store
schema 0.5 (0.1-0.4 still read); feature schema and production model
(`models/v5n`, SQ >= 5) unchanged.

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
