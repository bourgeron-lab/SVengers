# Changelog

## v2.1.4 — 2026-11-19

Bug fixes; model, features and feature schema unchanged. Calls from the CRAM
depth path (the paper's HG002 benchmark) are identical to v2.1.0.

- **MaxN no longer rejects polymorphic deletions masked by the panel of
  normals.** The MaxN hard filter fires when more than 70 % of a call's depth
  bins are masked, and that mask includes the bins the PoN masks as
  population-variable, which is exactly where common CNVs lie. With the PoN
  (the cohort fast path), 41 % of the true HG002 deletions of 1-10 kb (GIAB
  v5.0q) were tagged MaxN despite split and discordant reads. A deletion
  below 10 kb with junction evidence is now exempt from the masked-bin clause
  (the N-gap/exclusion clause still applies). Not above 10 kb, where the
  exemption passed megabase pericentromeric artefacts, and not for
  duplications, where it added 8,268 cohort sites (+62 %), mostly
  junction-only DUPs below 1 kb.
  - HG002 out of sample, cohort fast path: DEL F1 at 1-10 kb 0.465 -> 0.565
    (precision 0.966 -> 0.974); CMRG 0.455 -> 0.480.
  - 2,670 GHFC genomes: about 100 more PASS deletions of 1-10 kb per genome,
    the same in every wave stratum; in held-out trios 99 % of the added
    deletions are found in a parent.
  - 665 trios (`svengers bench cohort`): de novo DEL 1-10 kb median 6 -> 7
    per child, confirmed de novo 180 -> 182, Mendelian violation rate
    0.0234 -> 0.0236, duplications unchanged.
- **Cohort merge 10x faster, same output.** The clustering sweep compared
  every call with every open cluster of the last megabase; a cluster whose
  centroid ends before the sweep line can never take a member again and now
  leaves the scan at once, and cluster medians are kept on sorted lists
  instead of being recomputed after each member. Same clusters in the same
  order (checked on chr1, chr16 and chr21 of 2,670 genomes, and by a
  randomized test against the previous algorithm); the merge of 2,670
  genomes takes 2 minutes instead of 21.
- **Model card**: `docs/MODEL.md` gains the HG002 numbers of the cohort fast
  path, with and without the PoN.

## v2.1.3 — 2026-11-17

Documentation only; calls and scores are those of v2.1.0.

- **Locus-disjoint HG002 evaluation** (`docs/MODEL.md`) — the out-of-sample
  model trained without HG002 kept the HG003/HG004 labels, which are defined
  from HG002's truth set and false-positive loci. Two models trained without
  HG002 on half of the chromosomes each, HG002 scored on the chromosomes each
  never saw: deletion F1 0.515 against GIAB v5.0q (precision 0.817 instead of
  0.872), 0.843 at 1-10 kb, 0.717 at 10-100 kb, 0.552 against CMRG.
- **Model card corrected** — the GIAB genomes are NovaSeq PCR-free at about
  40x (not 30x), and the labels of HG003/HG004 are transmission labels
  against HG002 (`svengers label --parent-mode`), not their own truth sets.
- **Zenodo DOI** in the README and `CITATION.cff` (10.5281/zenodo.24308195,
  all versions).

## v2.1.2 — 2026-11-17

Citation metadata only; calls and scores are those of v2.1.0.

- **Archived on Zenodo** — releases are now archived with a DOI; the
  citation metadata (`CITATION.cff`) gains the author's ORCID and an
  abstract, which the archive uses as its description.

## v2.1.1 — 2026-11-17

Reproducibility outside GHFC; calls and scores are those of v2.1.0.

- **Container on ghcr.io** — every release is now published as
  `ghcr.io/bourgeron-lab/svengers:<version>` (and `:latest`), with the filter
  model (`/opt/svengers/models/v12b`) and the generic GRCh38 config
  (`/opt/svengers/configs/grch38.toml`) inside the image; the Docker and
  Apptainer recipes carry both.
- **The paper's reference files, rebuilt anywhere** — `configs/grch38.toml`
  points at the directory built by `workflows/build_grch38_resources.sh`,
  which downloads the GIAB GRCh38 FASTA, rebuilds the LCR mask
  (`selfchain-LCR.hg38.bed`, the slivar recipe: UCSC Self Chain with
  normScore > 90 plus LCR-hs38), takes the CNVnator exclusion list, the
  segmental duplications and `popCovStdev.bw` from the ClinSV reference data
  you download (its licence forbids redistribution), and checks every file
  against the checksum of the copy the paper used. Tested end to end from an
  empty directory: all six files identical. `docs/RESOURCES.md` lists
  sources, recipes, checksums, licences and what a missing track costs
  (without `popCovStdev.bw`, HG002 DEL F1 at 1-10 kb drops from 0.836 to
  0.801; without the LCR mask, overall DEL F1 drops from 0.515 to 0.503).
- **Documentation in English** — `docs/ROADMAP.md` and
  `docs/wave_qc_validation.md` translated, with the validation results of the
  frozen version (tool comparisons on HG002 and on 665 trios, confirmed de
  novo calls by affection status, provenance of the reference tracks);
  training and benchmark messages in English. The HTML report, family pages
  and the French text values of the exports are unchanged.

## v2.1.0 — 2026-11-17

**Method frozen for the paper.** From this version on, bug fixes only: no
change of model, filter or feature schema, so the paper's numbers hold.
Calls and scores are those of v2.0.2.

- **The filter model stays `models/v12b`, documented in `docs/MODEL.md`** —
  its training data (GIAB HG002/HG003/HG004, public; 16 children of family
  trios with wavy coverage, under controlled access), its labels and
  out-of-fold metrics. The same recipe trained on GIAB alone was evaluated
  as a fully public alternative: on HG002 out of sample it matches `v12b`
  (DEL F1 0.513 vs 0.515 against GIAB v5.0q; CMRG 0.536 vs 0.582), but on
  2,670 cohort genomes the wavy ones (wave score >= 0.20) get 3.3 times more
  deletions >= 10 kb just below PASS (3.2 to 10.4 per genome), the wave
  artefacts that cohort recurrence would rescue, and children lose about
  520 deletions of 1-10 kb that a parent carries. It is not shipped; the
  model card reports it as an ablation.

## v2.0.2 — 2026-11-16

- **scikit-learn pinned to 1.7** (`>=1.7,<1.8`, package and container
  environment) — the shipped filter models (`models/v12b` in production) are
  scikit-learn 1.7 pickles; a fresh install picked scikit-learn 1.9, which
  cannot read them (`No module named '_loss'`). Models trained with another
  scikit-learn version will need this pin moved with them.

## v2.0.1 — 2026-11-16

- **Models trained before 2.0 load on a fresh install** — their pickles name
  `vedette.filter.model.CalibratedModel`; v2.0.0 found that class only where
  the old `vedette` package was still importable, so `filter` with the
  production model (`models/v12b`) failed on a new install. While a model
  loads, `vedette.*` now resolves to the matching `svengers.*` module. Scores
  are unchanged (same class).

## v2.0.0 — 2026-11-16

Vedette is now **SVengers**, developed at https://github.com/bourgeron-lab/SVengers.
Results are unchanged: on a Phelan-McDermid cohort (319 genomes) re-filtered
from its v1.16.1 calls, every filtered VCF, the cohort VCF, the Mendelian
report and the family tables are identical but for the names.

- **New name everywhere it is written** — Python package and command
  `svengers`; per-sample raw calls `<S>.svengers.vcf.gz`; family exports
  `<FID>.svengers.{tsv,json,parquet,qc.tsv}`; cohort site ids
  `svengers_site_N`; VCF headers `##source=svengers-*` and
  `##svengers_coords`; export and model metadata `svengers_version`;
  Nextflow processes `SVENGERS_*` and parameters `--svengers_bin`,
  `--svengers_config`; container recipe `containers/svengers.def`.
- **What Vedette wrote stays readable** — the `vedette` command remains as an
  alias; raw calls named `<S>.vedette.vcf.gz` are found by `filter`,
  `regenotype`, `label` and the Nextflow `--calls_dir`; `##vedette_coords`
  and `source=vedette` headers keep their coordinate convention;
  `vedette_version` metadata, `<FID>.vedette.*` exports and the models
  trained so far load as before; `--vedette_bin` and `--vedette_config` are
  accepted; validations saved in the family pages under the old keys are
  still read. Evidence stores and the feature schema (hash
  f89b9e8890da996f) are untouched: no genome needs to be called again.
- **Nextflow** — `SVENGERS_BND` starts at 128 GB: at 64 GB it ran out of its
  6 h on 2,670 genomes, then took 44 min at 128 GB.

## v1.16.1 — 2026-11-16

The family tables and pages showed the wrong inheritance on merged rows.
Calling, filtering, genotyping and the Mendelian report are unchanged; only
`vedette report` changes.

- **Inheritance of a merged row** — a terminal or chromosomal event is cut
  into sub-sites, the confirmed de novo on the child's own call and the others
  only genotyped, and the row picks up small inherited polymorphisms. A vote
  by count over the sub-sites showed the 22q13.3 deletion of 16 children of a
  Phelan-McDermid cohort as `denovo_genotype_only`, a 10q26.3 deletion as
  biparental and a 17 Mb de novo duplication of chr18 as paternal (six
  paternal sites of 60-260 bp against three de novo of 2.7-7.9 Mb). The
  GQ-confident de novo labels now vote as one, weighted by length; when they
  cover at least as much of the row as the leading other label, the row shows
  the strongest of them (`denovo` > `denovo_weak` > `denovo_genotype_only`).
  Otherwise the vote stays by count, and a low-GQ de novo (`denovo_lowconf`)
  counts only when nothing else does: it no longer outvotes a confident
  inherited label. The rows with a confirmed de novo of the child shown under
  another label go from 50 to 0 in the Phelan-McDermid cohort (94 trios; all
  33 de novo events >= 1 Mb on a PASS call of the child now `denovo`, against
  3) and from 32 to 2 in dataset_paper (665 trios; the two left are the
  pseudo-autosomal rows of an XXY child, labelled by its karyotype). On
  dataset_paper 2,961 of 1.77 million family rows change label, mostly a
  low-GQ de novo replaced by the confident label beside it, and 964 "de novo
  < 1 kb" alerts go.

## v1.16.0 — 2026-11-16

Depth-only deletions outside repeats keep their PASS, and the cohort VCF flags
the genotypes held by split reads alone, which the trio tools now read as
missing. Model, store schema and feature schema unchanged; no retraining.

Batch is now measurable: each genome's flowcell, read from its read names,
gives trios whose parent and child sit on different flowcells, where a batch
artefact cannot pass for transmission (352 of 665 in dataset_paper).

- **`UNIQUE_DEL_DEPTH` override** — a depth-only deletion >= 10 kb (no
  junction evidence) outside repeats (segmental duplication and low-complexity
  overlap < 0.2), at CN <= 1, depth ratio <= 0.65, >= 90 % consistent bins,
  < 50 % masked bins, in a genome of wave score < 0.20, is PASS whatever its
  SQ (`INFO/SQ_OVERRIDE=UNIQUE_DEL_DEPTH`). In parent/child pairs on different
  flowcells, the parents' calls of this class at SQ < 5 are transmitted at 0.45
  against 0.07 in unrelated genomes of other flowcells (~85 % real, no worse
  than those already PASS). In repeats the class is ~60 % real at any SQ, which
  no threshold and no retraining on these features can sort: left to the
  model. The rule is off in a genome where it would fire more than 20 times
  (`unique_del_rescue_max_per_genome`): such a burst (up to 894 in one genome,
  chr4-heavy, against 1-11 elsewhere) is a genome-wide noise the wave score
  misses, and two such genomes carried nearly all the confirmed de novo a first
  version added. On the two production cohorts (2,989 genomes) it adds 1-2
  PASS deletions per genome in median (PASS DEL >= 10 kb +1.3 to +1.8 per
  genome in every wave stratum); on HG002 it adds 4, two of them true
  deletions of the GIAB v5.0q benchmark and two outside it. New `filter`
  settings `unique_del_rescue_size`, `_max_ratio`, `_min_consistent`,
  `_max_repeat`, `_max_mask`, `_max_wave`, `_max_per_genome`.
- **`FORMAT/FT=WeakSupport` in the cohort VCF** — a carrier genotype >= 1 kb
  with no discordant pair and a copy ratio that does not support its type (DEL
  RD > 0.75, > 0.62 from 10 kb; DUP RD < 1.25). The de novo left in wavy or
  batch-affected children are such genotypes (DEL 1-10 kb de novo: PE 0, RD
  ~0.75; inherited: PE ~10, RD ~0.6), up to thousands per child on one
  flowcell. `GT` is unchanged and the site filters still use it (the PASS site
  set is the same); `ped-annotate`, `bench cohort` and the trio benchmarks read
  a `WeakSupport` genotype as missing. On dataset_paper (665 trios) the median
  de novo per child goes from 19 to 15 (103 to 65 at wave >= 0.20), recurrent
  de novo sites from 1,401 to 1,012, for 2-8 % of the inherited genotypes
  >= 1 kb; the confirmed de novo go from 202 to 180 (26 split-read-only
  deletions of wavy children and GC-rich regions removed, 4 added). In the
  Phelan-McDermid cohort (94 trios) every child with the 22q13.3 deletion
  stays confirmed (135 to 129: redundant sites carrying another child's
  breakpoint). Mendelian violation and transmission rates unchanged. Thresholds
  in `cohort.metrics.TRIO_DEFAULTS` (`weak_*`).

## v1.15.1 — 2026-11-15

A clean terminal deletion keeps its PASS. Model, store schema and feature
schema unchanged; only the filter stage changes.

- **`TERMINAL_DEL_DEPTH` override** — a depth-only deletion >= 100 kb that
  reaches a chromosome end (within 300 kb of it, with <= 10 kb of callable
  sequence left beyond, as in `vedette bnd`), at CN <= 1, depth ratio <= 0.65,
  >= 90 % consistent bins and < 50 % masked bins, is PASS whatever its SQ
  (`INFO/SQ_OVERRIDE=TERMINAL_DEL_DEPTH`). A terminal deletion has no far
  junction, and model v12b, which weighs paired-end and junction evidence
  above all, scored such deletions at SQ ~1: the 134 kb 22q13.3 deletion
  (SHANK3) of eleven children of a Phelan-McDermid cohort stayed `LowSQ`, so
  their de novo was only "weak" and `vedette bnd` did not list the terminal
  loss. They are now confirmed de novo (DNQ 99). In trios, the parents'
  depth-only terminal deletions >= 100 kb are transmitted at 0.56-0.59 against
  0.10-0.20 for an unrelated child, whatever their SQ. On the two production
  cohorts (2,989 genomes) the rule touches 14 + 19 calls; the pre-gate and the
  de novo gate (by wave stratum) are unchanged, HG002 calls are identical, and
  no labelled training call is affected. New `filter` settings
  `terminal_del_rescue_size`, `_max_ratio`, `_min_consistent`, `_max_mask`.

## v1.15.0 — 2026-11-13

Model v12b replaces v5n, and the trios now confirm de novo calls and flag
sites: a confirmed de novo needs a PASS call of the child, no parental call
and a private site (lifted at a ClinGen dosage-sensitive locus), and the
cohort VCF gains the LowTransmission and DeNovoExcess flags. Sex-chromosome and
whole-autosome baselines come from the measured karyotype. Store schema and
feature schema unchanged; calling unchanged on euploid genomes.

- **Sex-chromosome baseline from the measured karyotype** — with `--pon`,
  `vedette call` sets the chrX/chrY ploidy from the karyotype (copy number on
  the panel's sex-informative bins, X heterozygosity, BAF) instead of the
  chrX < 0.75 / chrY > 0.25 threshold, which read an XXY as a man and an XXX
  as a woman. On dataset_paper, an XXY went from 168 PASS DEL / 4 DUP on chrX
  to 143 / 45 and an XXX from 109 / 133 to 129 / 45, both now like an XX
  woman (173 / 42); an XYY gets chrY 2. The PoN level is that of the profile
  used (0.5 male chrX/chrY, 1.0 female chrX and autosomes), and genomes with
  a sex-chromosome aneuploidy no longer enter the panel's sex profiles.
- **Whole-autosome aneuploidies get their own baseline** — a chromosome at 1
  or 3 copies over >= 85 % of its 1 Mb windows and confirmed by the BAF is
  called against that copy number. A trisomy 12 went from 128 to 13 PASS DUP
  on chr12 (the fragments of the trisomy) and 93 to 118 PASS DEL (deletions
  masked by the diploid baseline); a trisomy 21 lost its 15 DUP fragments of
  0.3 to 11 Mb (none left above 10 kb) and went from 107 to 130 PASS DEL
  on chr21. Checked at calling time on 87 genomes: the known trisomies 12
  and 21 only.
- Store attributes `karyotype`, `karyotype_note`, `aneuploidies`; sample VCF
  headers `##vedette_karyotype` and `##vedette_aneuploidy`. Re-call genomes
  with a sex-chromosome or autosomal aneuploidy to benefit; euploid genomes
  are unchanged (autosome calls identical, FILTER as before).

- **`vedette bench cohort`: the de novo gate before a model change** — compares
  cohort runs (the first is the baseline) on de novo calls per child at PASS
  sites with GQ >= 20 in the child and both parents, by type and size, and on
  recurrent de novo sites (de novo in >= 5 children). A median rise of 25 %
  and 2 calls per child in any class, or of 25 % and 20 recurrent sites,
  raises an alarm (exit code 1). On dataset_paper (665 children), the current
  code with model v5n is OK against production (median 27 -> 26), and a
  candidate model trained with IGV verdicts at weight 20 raises the alarm (DUP
  1-10 kb median 1 -> 9, DUP >= 10 kb 0 -> 3, recurrent sites 2,480 -> 3,994)
  although its per-genome PASS median on a 186-genome panel fell (DUP -9.7 %)
  and its extra child calls mostly overlapped a parent's call.
  The gate also runs within the children's wave strata (< 0.128, 0.128-0.20,
  >= 0.20, from the stores; >= 10 children): in the 27 children at wave >=
  0.20 that candidate went from 18 to 259 de novo DUPs each, and a second
  candidate (IGV deletion verdicts only), flat on every cohort median, from
  18 to 34 de novo DEL >= 10 kb.

- **`vedette bench pregate`: the cohort gate's quick first step** — on
  filter-only runs (`samples/<S>/<S>.filtered.vcf.gz`), compares per genome
  the carrier calls a candidate PASSes and its rescuable LowSQ calls (SQ >=
  the merge's `--rescue-min-sq`, which can rescue a site with no PASS
  member), by type and size, over all genomes and within the wave strata. A
  rise of the mean per genome of 25 % (PASS) or 50 % (rescuable) and 0.5
  call raises an alarm (exit code 1). On dataset_paper (2,670 genomes,
  minutes instead of ~2.5 h) it gives the cohort gate's verdict for the
  seven candidates the gate had judged, among them one that took the DUP
  1-10 kb PASS calls from 812 to 5,482 and the rescuable ones from 45k to
  136k (children's de novo DUP 1-10 kb median 1 -> 7) while the panel showed
  its DUPs falling by 11 %.

- **Genotyping: the integer CN reads edge bins with their flanks** — for
  events of at least one bin, each edge bin is a mix of the event and the
  flanks' copy number in proportion to the part of the bin it covers; read
  on whole bins, a 1.7 kb homozygous deletion inside a heterozygous one came
  out CN 1, and a heterozygous deletion that size CN 2. GT and GQ unchanged.
  The filter model learns on CN: calls genotyped before this version keep
  the old CN (re-run `call` or `vedette regenotype` them before filtering
  with v12b).
- **Model v12b replaces v5n** — v5n's labels and recipe (wave_score pinned
  at scoring) trained on that CN. On HG002 out of sample, DEL F1 0.514 vs
  0.511 (GIAB v5.0q) and 0.582 vs 0.563 (CMRG). On 20 dataset_paper trios
  the CN alone, under v5n, dropped 165 PASS 1-5 kb deletions a parent also
  carries; retrained, on the 186-genome panel (PASS per genome DEL +7.3 %,
  DUP -17.1 %) the children gain 737 DEL 1-10 kb, 97.7 % with a parent call
  (158 lost, 99.4 %), and 5,004 DEL < 1 kb at 60.9 % (859 lost, 78.7 %): the
  small deletions are its weak side. Cohort gate on the 665 trios of
  dataset_paper quiet in every wave stratum: de novo median per child 27 as
  in production, recurrent de novo sites 2,483 vs 2,480, PASS sites DEL
  59,727 vs 59,598 and DUP 14,083 vs 14,825; the pre-gate is quiet too.
- **`vedette ped-annotate --calls-dir`: a de novo needs a call of the child**
  — the cohort genotyper re-genotypes every site in every sample, so a child
  can be a carrier at a site other samples nominated while the caller found
  nothing in that child. With the per-sample filtered VCFs, a de novo whose
  child holds no call of the same type at the site (reciprocal overlap >= 50
  %, any FILTER) becomes `denovo_genotype_only`, outside `n_denovo`; the
  report gains `child_call`, the summary `n_denovo_genotype_only`. On
  dataset_paper, 57,705 of the 63,934 de novo at PASS sites (GQ >= 20) had
  no call of the child: median 7 per child instead of 27. A blind image
  review of 40 found no de novo among the 26 without a call and 10 among 14
  with a private PASS call and no parental call. The Nextflow pipeline passes
  its `samples/`; `bench cohort` gates on both labels (comparable with older
  runs) and prints the confirmed ones; the HTML report's de novo tables keep
  `denovo` only.

- **`ped-annotate --calls-dir`: a confirmed de novo needs a PASS call of the
  child, no parental call and a private site** — a call of the child was not
  enough: on dataset_paper 90 % of the de novo that kept one held a LowSQ
  call, mostly on sites de novo in several unrelated children (a genotyping
  miss in the parents, not a new mutation). `denovo` now also requires no call
  of the same type in either parent (any FILTER) and at most
  max(`--max-other-carriers` 2, `--max-other-carrier-frac` 0.005 x pedigree
  samples) carriers outside the family; the other GQ-confident de novo with a
  call of the child are `denovo_weak`. New report columns `parent_call` and
  `other_carriers`, summary column `n_denovo_weak`; `bench cohort` counts
  `denovo_weak` in the raw burden. dataset_paper (665 trios): 202 confirmed de
  novo (DEL 173, DUP 29).
- **No carrier limit at a ClinGen dosage-sensitive locus** — a pathogenic de
  novo recurs across families by construction: in the PMS cohort the terminal
  22q13.3 deletion over SHANK3 has 17 to 74 carriers per site outside the
  family, and every one became `denovo_weak`. With `--clingen-genes` /
  `--clingen-regions` (passed by the Nextflow ped step), the limit is lifted
  for a gene of HI 3 (DEL) or TS 3 (DUP), or a score-3 region covered >= 20 %,
  while at most `--max-other-carriers` of the outside carriers are parents in
  the pedigree: the recurrence must come from children (the small deletions of
  dataset_paper in 3q29, 22q11.2 or BRSK2 have 8 to 92 parent carriers). New
  report column `dosage_locus`. PMS (94 trios): 29 -> 113 confirmed de novo,
  84 at SHANK3 in 36 of the 47 children with the 22q13.3 deletion; dataset_paper
  unchanged. The 11 others hold a LowSQ call of the deletion (a 134 kb
  terminal loss seen by depth only), a weakness of the model.
- **`cohort-filter`: LowTransmission and DeNovoExcess trio flags** — with a
  pedigree, an autosomal site is `LowTransmission` when >= 10 trios with one
  heterozygous parent and the other 0/0 transmit it at a rate < 0.2
  (`--trio-min-opp`, `--low-transmission`; trios with GQ >= 20, `--trio-min-gq`), and `DeNovoExcess` when
  >= 3 carrier children, half of them or more, have both parents 0/0
  (`--denovo-excess-n`, `--denovo-excess-frac`). A de novo child whose copy
  ratio supports the type (DEL RD <= 0.75, DUP RD >= 1.25) does not count, so
  true recurrent de novo deletions (22q13.3, 11q13.3) stay PASS. INFO
  `TR_OPP`, `TR_RATE`, `DN_TRIO`, `DN_FRAC`, `DN_NODEPTH`. On dataset_paper
  they remove 8,232 DEL and 581 DUP sites that were otherwise PASS; with the
  other v1.15 changes, PASS sites DEL 59,598 -> 51,517 and DUP 14,825 ->
  13,233, and de novo per child (raw burden) median 27 -> 19, 197 -> 103 in
  the wavy children (wave >= 0.20). Cohort gate quiet on both cohorts.

## v1.14.0 — 2026-11-07

`vedette bnd` calibrates standalone inversions on trios: same-strand
junctions pass at 12 pairs except in genomes with chimeric-library noise, and
common inversions are kept with FILTER `Common`. Store schema, feature schema
and model (v5n) unchanged; calling and genotyping unchanged.

- **`vedette bnd`: same-strand junctions pass at 12 pairs** — standalone
  inversion and fold-back junctions with split reads at both ends are PASS
  at the translocations' threshold (12 discordant pairs) instead of all
  filtered `Inversion`. On 637 dataset_paper and 94 PMS trios, 1.32 and 1.35
  per child are in a parent and 0.031 and 0.074 in the child alone; with an
  unrelated family's parents, 36 % match (population inversions). The ~1
  child-only junction per trio came from a few genomes: a median genome has
  7 kept same-strand junctions (PMS 11), 67 of 2,670 have more than 500, and
  their same-strand junctions are now filtered `NoisyGenome`.

- **`vedette bnd`: common inversions are kept** — a same-strand junction in
  1 % of the run or more was dropped as recurrent; with 12 pairs or more it
  is now kept (flag `common`, FILTER `Common`). On 633 dataset_paper trios,
  recurrent same-strand loci split by transmission: at 605, a child without
  a carrier parent carries the junction at most a quarter as often as the
  run (17q21.31 H2: 42 % of genomes, 30 pairs), at 942 as often as anyone
  (artefacts, a median of 4 pairs). Of the junctions with 12 pairs or more,
  74 % are at the first, 4 % at the second. ~110 per genome in the TSV
  (PMS ~140), ~60-70 in the VCF.

## v1.13.0 — 2026-11-07

`vedette bnd` reads translocations, inversions and the mechanism of terminal
losses from the evidence stores, and the report shows them. Cohort genotyping
reads the store's original BAF and ignores the hets of paralog-dense events.
De novo candidates are sorted by the child's genome quality. Same store
schema, feature schema and model (v5n).

- **Genotype likelihood: hets of paralog-dense events are ignored** — when
  an event holds more hets than a diploid stretch of its length yields
  (Poisson at twice the sample's contig het rate, tail probability
  below `[genotype] baf_paralog_p` = 0.05), its BAF term is dropped: those
  hets come from a paralog's reads whatever the copy number. Cohort
  genotyping of dataset_paper chr12 on the rewritten BAF: 3,653 de novo
  without the check, 2,921 with it (2,840 on the original BAF), Mendelian
  violations 280,264 -> 271,195 (268,078); on the original BAF the check
  is neutral (2,862 de novo, 267,883 violations). `vedette call` does not
  set the contig het rate, so calling-stage genotypes are unchanged;
  `vedette regenotype`, which reads the rewritten BAF, applies the check.

- **Report: a terminal loss a parent also carries is marked inherited** —
  "héritée de <parent>" when a parent has the same terminal loss (same
  contig, reciprocal overlap >= 0.5). With the callable-tail criterion, the
  only terminal losses with identical bounds in two genomes of PMS and
  dataset_paper are a mother's and her child's (chr5:0-141 kb, chr11:134.9
  Mb to qter), both PASS sites carried by those two genomes only.

- **Cohort genotyping reads the store's original BAF** — `vedette
  genotype` reads `<contig>/baf_track` (the GHFC track, QUAL >= 30) when
  `vedette baf` rewrote the store, unless `[genotype] rewritten_baf =
  true`; `vedette regenotype` still reads the rewritten BAF. The likelihood
  counts every het inside a deletion as a mismap (log 0.01 each, up to 50),
  and in paralogous loci the rewritten BAF has ~20 hets per kb: a common
  1.2 kb deletion at chr12:31.12 Mb, genotyped 0/1 in 98 % of dataset_paper
  with identical PE/SR/depth, fell to 40 % carriers, escaped CohortArtifact
  and HighCarrierFrac and gave de novo calls wherever the parents fell to
  0/0. Re-genotyped on the rewritten stores with v5n, dataset_paper had a
  median of 42 PASS de novo (GQ >= 20) per child against 27 in production,
  16,800 of them at 1,503 sites production flagged as artifacts. On chr12
  (2,670 genomes): 3,025 de novo in production, 3,653 with the rewritten
  BAF, 2,840 with the current code on the original BAF (de novo rows at
  sites recurrent in >= 5 children: 996, 1,669, 870).

- **Report: de novo candidates sorted by the child's genome quality** — a
  "qualité enfant" column names what weakens the child's genome (wavy
  coverage, old genome, excess of carried sites), and each family's table
  lists the de novo of good genomes first, then by impact: on the IGV
  review of sub-threshold rows (2026-10-13), calls were true 1 time in 10
  to 13 on old or wavy genomes, 15 in 28 elsewhere.

- **`vedette bnd`: translocations, inversions and the mechanism of terminal
  losses** — a run-level pass on the evidence stores (no CRAM). Inter-
  chromosomal pairs make translocation junctions and same-strand pairs of
  1 kb or more inversion junctions, reciprocal or fold-back (pairs on one
  side: one junction of an inverted duplication), each end placed on its
  split-read pileup (`<S>.bnd.tsv`; the raw scan is cached in
  `<S>.junctions.tsv`). Over the run, a locus joined to 10 distinct loci or
  more is a repeat (chr2:32.9 Mb: ~24,000 in PMS; a translocation's ends:
  1-3) and two loci joined in 5 genomes or more (or 1 %) an artefact or a
  common variant; counting genomes per locus does not work, a median locus
  has junctions in 7 PMS genomes. Each deletion of 100 kb or more reaching a
  contig end, with at most 10 kb of callable sequence past it, then gets its
  mechanism (`<S>.terminal.tsv`): translocation (with the partner's gain when
  the sample has it), inv-dup (a fold-back on the retained side), healed
  telomere (TTAGGG clips within 30 kb on the retained side, or a
  subtelomeric junction), a junction to a repeat within 2 kb of the
  breakpoint (rDNA, acrocentric short arms, decoys: sequence the reference
  does not place), or unknown. A loss called in abutting pieces breaks at its inner piece: 10
  PMS losses seemed to break at 49.32 Mb, the outer piece of losses from
  43.6-49.1 Mb. The report has a "Pertes terminales" section and flags the
  family rows of these losses; a parent with the same junction is a
  balanced carrier. Junctions with split reads at both ends are also written
  as VCF breakend pairs (`<S>.bnd.vcf.gz`, MATEID / EVENT / KIND); PASS
  needs 12 discordant pairs: on 100 dataset_paper trios, junctions in the
  child and in neither parent number ~180 per trio with 4-7 pairs, 2.9 with
  8-11 and 0.2 with 12 or more. Same-strand junctions stay `Inversion`
  (~1 per trio even at 12 pairs). On the 63 PMS terminal losses: 9
  translocations (6 with the partner's end gained), 5 inv-dup (3 with the
  adjacent duplication called), 29 healed telomeres, 6 junctions to a
  repeat, 14 unknown (5 of them mosaic).
  The scan takes ~20 s per genome (pileups made once per contig).

## v1.12.0 — 2026-11-06

Large depth segments are no longer cut by one or two outlier bins, `vedette
train` takes per-row weights, and the mosaic scan places its edges at 1 kb
and keeps a clone in one row. Same store schema, feature schema and model
(v5n); the segment bridging takes effect on the next `vedette call`.

- **`vedette call`: a large segment broken by a short interruption is
  joined back** — the Gaussian depth emission has light tails: one 1 kb bin
  at ratio ~1.3 inside a het deletion costs ~100 nats under CN1 against 37
  for two state switches, and the 2 kb CN0 of a common deletion on the
  remaining allele (chr22:49,319,854-49,321,519, carried by 84 of 319 PMS
  genomes) does the same. A 6.8 Mb terminal 22q13 loss came out as 8
  abutting DEL calls (2 of them LowSQ), and 10 PMS genomes had their
  terminal deletion cut at exactly 49,322,000. Two same-CN segments are now
  joined when the bins between them number at most
  `segment.bridge_max_gap_bins` (3) and at most `bridge_max_gap_frac`
  (0.05) of the shorter side, and their mean ratios agree within
  `bridge_max_ratio_diff` (0.15); the interrupting segments stay as their
  own nested calls, and the Viterbi path is unchanged. Replayed on the
  production VCFs: on PMS 410 joins (869 calls -> 410), 26 of them with a
  PASS piece, among which the 22q13 deletions of 14 genomes; on
  dataset_paper 4,868 joins (10,627 -> 4,868), 127 with a PASS piece,
  mostly aneuploidies (trisomy 12, sex chromosomes) cut in dozens of pieces.
  Small events are untouched (a 1-bin gap needs 20 bins on each side).
  Takes effect on the next `vedette call`; store and feature schemas and
  the model are unchanged.

- **`vedette train`: optional `weight` column in the labeled TSVs** — the
  number of calls a row stands for (default 1), used by the trees and the
  isotonic calibration; the manifest's metrics still count rows and record
  `weighted_rows` / `total_weight`. A few curated verdicts (129 IGV labels
  on depth-only calls) can then weigh against ~50,000 truth and
  transmission labels.
- **Mosaic scan: inner edges refined to the store bin (1 kb)** — segments
  were cut at 1 Mb windows. Each inner edge is now placed where the 1 kb
  depth steps from the next window out to the segment's level, within a
  window of the old edge. A PMS proband's 22q13 loss at f = 0.18 moves from
  42.0 Mb to 41,993,000, 6 kb from its telomere-healing junction (one split
  read); a 38 % loss seen in IGV at 47,730,170 comes back at 47,741,000.
  The edge stays where it was when there is no step to place: the next
  window out is masked (a centromere) or within a quarter of the segment's
  shift from its level (a gain running on past the window), or the step
  lands on a limit of the search. On PMS and dataset_paper the same 20 and
  26 segments, edges moved by a median of 0.28 Mb (at most 0.72 Mb).

- **Mosaic scan: segments no longer cut at their straddling windows** —
  successive 1 Mb windows differ by ~0.01 of the diploid level, so every
  step is significant and the window straddling a breakpoint, kept with two
  neighbours to make the 3-window minimum, formed a blended segment the
  homogeneity test dropped. Such windows are now left out of the
  segmentation (the 1 kb refinement finds the breakpoint inside them), a
  blended segment's end windows go back to the neighbour whose level they
  share, and adjacent segments within 0.1 copies are joined. On PMS and
  dataset_paper: no segment lost or gained; a 22q13 duplication goes from
  44-48 Mb to 41.5-49.3 Mb, and one clone or trisomy split into 2-3 rows
  is one row (a clonal 17q gain, chromosomes 8, 12 and 21): 57 rows -> 52.

## v1.11.0 — 2026-11-05

- **Mosaic scan: a wavy genome needs a fraction of at least 0.3** — depth
  waves of ±15 % read as 15-25 % clones. On the BAF rewritten by `vedette
  baf` (twice the heterozygous sites, a tighter window noise), 9 segments at
  f = 0.15-0.26 passed the BAF test on genomes with a wave score (without
  the segment's contig) of 0.134-0.210, mostly on GC-rich chromosomes. With
  the floor, the 2,960 rewritten genomes give the 22 segments of the 9
  carriers found on the old BAF, plus one 22q13 loss at f = 0.18 on a
  genome at 0.10, the BAF exactly as expected. On the rewritten BAF the
  known mosaics' BAF excess matches its expectation (0.146 for 0.146,
  against 0.108 for 0.157 before).

## v1.10.0 — 2026-11-05

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
