# SVengers roadmap

State as of 2026-09-08: v0.1 core working end to end (call → filter →
merge → genotype → cohort-filter → ped-annotate), validated on the `test`
trio chr20 (breakpoints identical to ClinSV, MVR 0.6%), 34 unit tests.

## Next steps (by priority)

### State as of 2026-11-17 (v2.1.0): method frozen for the paper, then validation

Vedette has been called **SVengers** since v2.0.0 (repository
bourgeron-lab/SVengers; `vedette` remains an alias of the command; what
Vedette wrote remains readable, see `legacy.py`). v2.0.2 is in production on
dataset_paper (2,670 genomes) and PMS (319): refiltering identical to v1.16.1
apart from names, pre-gate and cohort gate calm, CI green. Since the state of
2026-11-10: v1.15.1 (clear terminal loss kept PASS), v1.16.0
(`UNIQUE_DEL_DEPTH`: depth DEL in unique sequence kept PASS, rule switched
off in a genome where it would fire more than 20 times; FORMAT/FT
`WeakSupport`), v1.16.1 (inheritance label of merged family rows),
v2.0.1-2.0.2 (models pickled before 2.0 readable on a fresh install,
scikit-learn pinned to 1.7).

**Method frozen for the paper: v2.1.0 (2026-11-17).** Since the freeze:
bug fixes only; no change to the model, filters or feature schema,
otherwise the paper's numbers move.

**0. Paper model: v12b (decided on 2026-11-17).** v12b is trained partly
on 16 GHFC children (71% of the labels, all of the wavy genomes). Its
weights are published with the code, the training data are declared as
controlled access (`docs/MODEL.md`). Retraining on public data only
(**v14p**, GIAB HG002/3/4, same recipe, same features and labels;
`runs/_v14p/`) was evaluated and not retained:
- HG002 out of sample: DEL F1 0.513 versus 0.515 (GIAB v5.0q), 1-10 kb
  0.843 versus 0.836, CMRG 0.536 versus 0.582;
- panel of 186 genomes: PASS per genome DEL −2.8%, DUP −25% (mostly
  < 1 kb); in the trios, children lose ~520 DEL of 1-10 kb, 97.9% of them
  carried by a parent;
- **dataset_paper pre-gate in ALARM**: wave stratum ≥ 0.20 (62 genomes),
  rescuable DEL ≥ 10 kb 3.16 → 10.44 per genome; cohort gate not run.
Without wavy genomes in training, the model no longer rejects coverage
waves. v14p becomes the paper's ablation. It remains to be confirmed on the
governance side that releasing the weights is compatible with the cohorts'
access agreements. Only remaining 100% public route: public wavy trios
(1000 Genomes high coverage, NYGC), ~1 TB to download, outside the freeze.

**Validation on the frozen version, by cost/value ratio.**
1. **Orthogonal sensitivity**: the 22q13 (SHANK3) deletions of PMS
   patients known from another technique (array, ACPA, MLPA), with their
   sizes; data to be provided by Alexandre. Gives a real sensitivity by
   size and the accuracy of the breakpoints.
2. **Comparisons without re-running anything: done (2026-11-17)**,
   `runs/_cmp/`. Still remaining: the depth callers on HG002 (CNVnator,
   GATK-gCNV), only if the full paper is chosen.
   - **HG002** (`runs/_cmp/hg002/`, 2026-11-10/12; v12b trained without
     HG002, PASS; ClinSV run on HG002) — F1:

     | | SVengers | Manta | dysgu | DELLY | ClinSV |
     |---|---|---|---|---|---|
     | DEL all sizes (v5.0q) | 0.514 | 0.536 | 0.541 | 0.423 | 0.395 |
     | DEL 50-300 bp | 0.396 | 0.458 | 0.454 | 0.259 | 0.191 |
     | DEL 1-10 kb | 0.831 | 0.769 | 0.778 | 0.800 | 0.691 |
     | DEL 10-100 kb | 0.718 | 0.595 | 0.500 | 0.456 | 0.534 |
     | Tandem DUP (DUP→INS) | 0.145 | 0.057 | 0.004 | 0.032 | 0.020 |
     | DEL CMRG | 0.582 | 0.652 | 0.582 | 0.400 | 0.426 |

     SVengers ranks first on DEL from 1 to 100 kb (precision 0.963 at
     1-10 kb, 1.000 at 10-100 kb) and last among the junction callers below
     300 bp.
   - **dataset_paper trios** (`runs/_cmp/trios_tools.py` →
     `trios_tools.tsv`, 665 children, individual PASS calls of each tool,
     autosomes): for each PASS call of the child, a call by the same tool
     in a parent (RO ≥ 0.5, any FILTER; also PASS only); "unexplained" =
     without a parental call (de novo candidates, almost all false);
     transmission of PASS calls private to one parent (expected ~0.5).

     | DEL | per child | with parent | unexplained/child | transmission |
     |---|---|---|---|---|
     | SVengers 1-10 kb | 276 | 99.3% | 2 | 0.55 |
     | ClinSV 1-10 kb | 878 | 95.3% | 40 | 0.53 |
     | SVengers 10-100 kb | 18 | 98.4% | 0.3 (mean) | 0.51 |
     | ClinSV 10-100 kb | 564 | 78.9% | 122 (mean) | 0.40 |
     | WisecondorX 10-100 kb | 5 | 59.7% | 3.5 (mean) | 0.41 |
     | SVengers ≥ 100 kb | 1 | 80% | 0.2 (mean) | 0.35 |
     | WisecondorX ≥ 100 kb | 1 | 24% | 15.7 (mean) | 0.13 |

     Genomes with wave ≥ 0.20 (27 children): SVengers unexplained DEL
     1-10 kb 3 per child (2 elsewhere); ClinSV DEL 10-100 kb 331 on
     average; WisecondorX DEL ≥ 100 kb 249.
   - **To present as is**: DEL < 1 kb, SVengers 392 unexplained per child
     versus 178 for ClinSV; DUP ≥ 1 kb, almost none PASS at the individual
     level (calls present, SQ below the threshold; the cohort stage
     recovers them by recurrence), while ClinSV outputs hundreds per child;
     with the parent in PASS only, SVengers DEL 1-10 kb 87.7% with parent
     versus 92.1% for ClinSV (the SVengers parental call is often LowSQ).
3. **Out of sample**: HG002 with a model trained without HG002;
   Platinum Pedigree (plan in "Benchmark harness"; download with
   Alexandre's approval, sizes announced beforehand).
4. **Reproducibility**: generic config, GRCh38 resource bundle, public
   container (already built in CI), Zenodo DOI; PyPI and bioconda only
   with Alexandre's approval. Reference tracks (2026-11-17, `runs/_repro/`):
   - `selfchain-LCR.hg38.bed` rebuilt byte-identical (md5
     b4406771c8fe431d114430280b39bb60) with the slivar wiki recipe
     ("finding deletions in parent-child duos", hg38): UCSC hg38
     `chainSelf` (2020-06-21) with `normScore > 90`, target side and raw
     query coordinates (no minus-strand conversion, so ~41 Mb of the mask
     sit at mirrored positions; kept as is, v12b learned on it) + lh3
     `LCR-hs38`, `chr` stripped, `bedtools merge`
     (`runs/_repro/lcr_slivar.sbatch`). Public sources: ship the recipe and
     the expected md5, or the file itself with its citations.
   - `popCovStdev.bw` comes from ClinSV refdata-b38, whose licence forbids
     redistribution: document the download (39 GB archive, research use).
   - Cost of a missing track (v12b without HG002 on HG002; 186-genome
     panel): without popCovStdev, DEL 1-10 kb F1 0.836 -> 0.801, PASS DEL
     +12% (noise < 1 kb), DUP -27%, ~80 real 1-10 kb deletions lost per
     child; without the LCR mask, DEL F1 0.515 -> 0.503, PASS DEL -12%,
     DUP +26%. A substitute track is not an option during the freeze.
   - Still to check: the segmental-duplication and exclusion BEDs are
     copies from ClinSV refdata; their public origins (UCSC
     genomicSuperDups, speedseq) are not yet verified byte-identical.
5. **Writing**: application note or full paper, to be decided (sets the
   scope of 2 and 3). All de novo numbers are to be recomputed on the
   frozen version at the "confirmed de novo" level (PASS call of the child,
   no parental call, private site; already measured by `bench cohort`).

**Outside the freeze** (after the paper, or in parallel without touching
production): fragmented calls of large deletions and complex events (`bnd`,
points 3 and 5 below); DUP weaker than DEL (to present as is in the
paper); recurrence rescue of 1-10 kb DUP; CNN second opinion
(1sexdecies); SVD on the panel of normals; **modern, redistributable
reference tracks for the next model**: replace `selfchain-LCR.hg38.bed`
(2016-era slivar recipe, built for `duo-del`, ~12% of the mask at
mirrored positions) by the GIAB GRCh38 stratifications v3.x (tandem
repeats and homopolymers, low mappability, segmental duplications,
all difficult regions; public domain, versioned, the same regions the
benchmarks stratify on), consider an SV-oriented LCR mask (e.g. the
2025 GRCh38 set of 111,067 regions, 35.4 Mb, from "Challenges in
structural variant calling in low-complexity regions") and the ENCODE
blacklist v2, and replace `popCovStdev.bw` by an open population-depth
track (e.g. built from public high-coverage genomes). Then retrain and
run the full deployment protocol (HG002/CMRG, panel, trios, pre-gate,
cohort gate by wave stratum); the goal is a model with no
non-redistributable resource.

### State as of 2026-11-10 (v1.15, unreleased) and priorities

On `main` since v1.14.0: `bench cohort` (de novo gate, by wave strata),
`bench pregate`, ploidy from the measured karyotype (sex chromosomes and
whole autosomes), integer CN read at the edges with the flanks, model
**v12b** by default, `ped-annotate --calls-dir` (`denovo_genotype_only`).
Production re-call in progress (dataset_paper then PMS): CN injected from a
re-genotyping, full re-call of the genomes whose karyotype changes, cohort
run in a separate directory, pre-gate and gate against production, then
swap (archives `_avant_v115/`) and refresh of the exports.

**1. Cohort de novo calls were mostly forced genotyping.** The cohort
genotyper re-genotypes every site in every sample: a child can be a carrier
at a site proposed by others without the caller having seen anything in
that child. dataset_paper (665 trios): 57,705 of the 63,934 de novo (PASS
sites, GQ ≥ 20) without any call of the child at the site, all filters
combined. Blind review of 40 cases (store images + reads):
- PASS call of the child, no parental call, private site (14): 10 true
  de novo (including 2 fragments of a larger event), 2 inherited CNVs
  (probable mosaic parent, carrier mother on chrX), 1 uncertain, 1 false;
- without a call of the child, private site (13): no de novo at the site,
  3 fragments of a larger event covered by the child's call;
- without a call of the child, recurrent site (13): 13 false.
False ones: nothing visible (the majority), GC coverage waves, multicopy
pericentromeres, HLA-DRB, an XXY read as XY, a small common homozygous DEL
wrongly sized at the site. Estimate: ~0.14 true de novo per child, versus
27 "de novo" per child displayed until now. Done: `denovo_genotype_only`
(median 7 per child remaining `denovo`, a necessary condition, not a
sufficient one: inherited variants mis-genotyped in a parent, recurrent
sites).
Remaining:
- a "confirmed de novo" level: PASS call of the child, no parental call at
  any filter, private or nearly private site; to be validated on a second,
  broader review;
- sites rescued by recurrence (RESCUE=RECURRENT) carry 83% of the de novo
  1-10 kb DUP: a flag or separate handling;
- the `bench cohort` gate therefore measures an artifact burden: useful for
  relative comparisons between models, it is not a de novo rate;
- every de novo number intended for the paper is to be recomputed.

**2. Models after v5n: what was rejected and why.**
- v12a (IGV ×20), v12a5 (×5), v12ad (IGV DEL only): de novo DUP 1-10 kb
  median 1 → 8-9 per child (v12a, v12a5); v12ad neutral in median but
  DEL ≥ 10 kb of the 27 children with wave ≥ 0.20 from 18 to 34.
- v13a-d (cleaned transmission labels of the 16 wavy families):
  transmission labels reward systematic artifacts, which are "transmitted"
  like real variants. v13c (< 1 kb labels of the wavy families removed):
  best HG002 F1 (0.520) but de novo DUP 1-10 kb 1 → 7; cause: "split-read
  only" DUP (junction, PE 0, CN 2, ~1.5 kb) labeled true 262 times / false
  49 in these families, versus 7 / 90 in GIAB, and SQ raised for LowSQ
  calls (recurrence rescue × 3.7). v13d (< 1 kb negatives put back): even
  worse in wavy children (DUP ≥ 10 kb 1 → 33). Lineage abandoned.
- v12b retained: v5n recipe on the new CN. HG002 out of sample F1
  DEL 0.514 (v5q) / 0.582 (CMRG) versus 0.511 / 0.563; panel DEL +7.3%,
  DUP −17.1%; gate calm in all strata. Weak point: DEL < 1 kb gained
  5,004 at 61% with a parental call, lost 859 at 79%.
- Lesson: the panel of 186 genomes said the opposite of the cohort for DUP
  three times (DUP per genome down, de novo DUP up); it has too few
  complete trios.

**3. Model deployment protocol (now mandatory).**
HG002 out of sample and CMRG; panel PASS per genome; trios "child call with
parental call"; then a `--steps filter` run of the whole cohort and
`bench pregate` (a few minutes; same verdict as the gate on the 7
candidates assessed); only then the full run and `bench cohort` by wave
strata.

**4. Publication plan: status** (detailed plan below, 2026-10-14).
1. Out-of-sample evaluation: partly (v12b numbers without HG002); a truth
   set independent of GIAB is missing (Platinum Pedigree, HGSVC).
2. Comparisons with other tools: partly. HG002 GIAB v5.0q with Manta,
   DELLY, dysgu (2026-09-09, `resources/giab/comparison_final`, to be
   re-scored against v12b); ClinSV (lumpy + CNVnator) and WisecondorX
   exist for ~4,000 GHFC genomes including dataset_paper (trio comparison
   possible without re-running anything). Missing: depth callers on HG002
   (CNVnator, ClinSV; GATK-gCNV, Canvas not installed), trio metrics
   compared across tools, Platinum Pedigree.
3. Orthogonal validation (array, MLPA, PMS SHANK3 deletions): not done,
   only IGV reviews.
4. Coverage waves: median PoN in production, SVD set aside, `wave_score`
   still neutralized in the model.
5. Scope: DUP clearly worse than DEL; overall DEL recall 0.36 on HG002
   (DEL < 300 bp, segmental duplications); de novo to be redone (point 1).
6. Reproducibility outside GHFC: not done (model partly trained on 16
   private GHFC children, GHFC config, no resource bundle, no container, no
   DOI); a blocker for publishers.
7. Writing: not started.

**5. Notes.**
- PAR1/PAR2 of males at CN 1 in the store: not a bug. The CRAMs have the
  chrY PAR masked (no reads), the raw depth is indeed at 2 copies, and the
  PoN normalization (male profile) brings all of chrX to the one-copy
  level, consistent with the genotyping. The outliers are real anomalies:
  XXY, XYY, XXX, mosaic loss of Y in a father, familial PAR1 gain.
- The image prototype `runs/_plot/plot2.py` (untracked) compares the
  PoN-normalized depth to the raw panel profile: false deviations in male
  PAR and at polymorphic loci. To be fixed before the next review.
- Tests: on a compute node with `OMP_NUM_THREADS=1` (~1 min), then the 4
  Truvari tests on the login node (SIGILL on the compute nodes).

### Translocations, inversions and terminal losses (`svengers bnd`, 2026-11-06)

Done (branch `bnd-terminal`): inter-chromosomal and same-strand junctions
read from the stores (no CRAM), cohort flags by distinct partner loci
(≥ 10: repeat) and by recurrent locus pair (≥ 5 genomes or 1%), mechanism
of each terminal loss (translocation, inv-dup, neo-telomere, junction to a
repeat, undetermined), "Terminal losses" section of the report with the
balanced carrier parent, breakend VCF (PASS: inter-chromosomal ≥ 12
pairs). PMS: 63 terminal losses, 9 translocations (6 with the partner
gain), 5 inv-dup, 29 telomeres, 6 repeats, 14 undetermined (including 5
mosaics). 100 dataset_paper trios: child-only junctions 0.13 per trio at
the PASS threshold.

Still to do:
1. ~~Nextflow workflow step~~: `SVENGERS_BND` before the report (tested
   on the public trio, 2026-11-06).
2. ~~Isolated inversions~~ (DONE): the ~1 per trio came from a few
   genomes with a chimeric library. Same-strand junctions placed at ≥ 12
   pairs, 637 dataset_paper trios / 94 PMS: 1.32 / 1.35 per child in a
   parent, 0.031 / 0.074 in the child only (parents from another family:
   36%, population inversions). PASS at the same threshold as the
   translocations, except in a genome with > 500 same-strand junctions
   (`NoisyGenome`, 67 / 2,670). Common inversions (≥ 1% of the run): a
   recurrent same-strand junction at ≥ 12 pairs is kept (`common`, FILTER
   `Common`); over 633 trios, 74% of these junctions are at loci that a
   child without a carrier parent carries at most 4 times less often than
   the run (17q21.31 H2: 42% of genomes, 30 pairs), 4% at loci carried at
   random (artifacts, ~4 pairs). ~110 per genome (PMS ~140). HG002
   benchmark to do (inversion truth).
3. Complex events: chain the junctions of a same genome (inv-dup without a
   called DUP, three-partner translocations, chromothripsis), and link
   junctions to intra-chromosomal DEL/DUP breakpoints.
4. The 14 undetermined losses in PMS: imprecise depth breakpoints (~20 kb;
   strong clip pileups there are often 10-20 kb away), five mosaics with
   few reads.
5. Fragmented calls: a large deletion comes out as adjacent pieces (last
   piece at 49.322 Mb in 10 PMS genomes); `bnd` chains them, the calling
   step needs to be fixed (investigation in progress).

### State as of 2026-11-05 (v1.7.0) and priorities

In production: v1.7.0, `models/v5n`, SQ >= 5, median PoN, sex chromosomes
handled end to end, impact read from the ClinGen dosage map (v1.5.0). The
production cohort VCFs keep the old coordinate convention until the merge is
re-run.
Done from 2026-11-02 to 2026-11-04 (block 1 below):
"events" pages for a sex-chromosome aneuploidy ("XXY karyotype"),
`ChrYNoDepth` extended to chrY sites < 1 kb (55 chrY PASS sites on
dataset_paper, father/son discordance of DELs ≥ 1 kb 10%, the real 6.9 kb
DUP transmitted 16/16 keeps its PASS), GC safeguard for the karyotype (two
false X/Y "mosaics" on very wavy genomes, GC contrast 1.06-1.23, moved to
"?"), old male stores in the archives deleted (862 GB). Outputs exported to
the GHFC tree (`samples|families|cohorts/…/svs/vedette/`, untracked scripts
under `runs/_export/`): to refresh after every production re-run.
Pedigrees corrected from the genomes (karyotype, relatedness on the BAF
sites) next to the originals, `*.pedigree.fix.tsv` + log: PMS 33
corrections (77 -> 94 trios analyzed, re-run), dataset_paper 4 corrections
(no trio changes, to include in the next production re-run).

genetics-viz issues (SVengers tab of the GHFC curation tool, #5-#11):
v1.4.0 DONE — VCF coordinates to the standard (`POS` = padding base, `END` =
last base; before: both one base to the right, and the merge added +1 to
the site), machine export of the family pages (`<FID>.svengers.tsv/.qc.tsv/
.json`, docs/EXPORT.md: #5, #7, #10). v1.7.0 DONE — non-PASS sites carried
by <= 20% of the cohort in the export, not on the page (#10, still open; the
others, on dataset_paper, are 88% carried by more than half of the
cohort). Remaining, to do together during the BAF recomputations
(point 3): DONE (v1.9.0) — `svengers labels --from-genetics-viz` (#11;
`validations/vedette.tsv` absent to date), and in a single change of the
export format: each member's own call with SOURCE, precision,
SQ/PE/SR/CG/RD/CNQ and the site INFO (#8), `site_key` + `svengers site-map`
(same breakpoints, then mutual best pairs at reciprocal overlap
>= 0.5, then split or merged sites: #9) and long Parquet (#6).

Priorities:
1. **Independent labels for depth** — IGV batch v5 REVIEWED (2026-11-05):
   141 depth-only calls (children of dataset_paper trios),
   stratified by type × size × score zone × apparent transmission,
   blinded (`runs/validation/`, untracked; analysis `analyse_v5.py`).
   84 true, 47 false, 10 uncertain. PASS: 30/30 decided true. Below the
   threshold: 54 true / 47 false, whatever the zone (grey 55%, low 53%);
   SQ no longer ranks anything there (AUC 0.44 DEL, 0.54 DUP). What
   separates them: apparent transmission (71-75% true versus 31-33% for
   apparent de novo; the review saw the parents' profiles), wavy genome
   (33% versus ~62%), DUP size (80% ≥ 100 kb, ~40% below).
   **v5n loses real large DUPs**: 12/15 true ones below the threshold (8/8
   inherited), breakpoints too short on the 3 corrected ones. Follow-ups:
   a. impact score with the ClinGen dosage map — DONE (v1.5.0):
      a de novo dup15q (CN 4, affected child) was at +0.15;
   b. rescue of large DUPs below the threshold — DONE (v1.6.0): the "SOUS-SEUIL GRANDE
      DUP" report label, LowSQ DUPs at CN ≥ 3, GQ ≥ 90, joined within
      the same carrier (gap ≤ 10% of the span, ≤ 1 Mb), ≥ 100 kb, private
      to the family, excluding wavy/IP genomes and X/Y of an atypical
      karyotype; the dup15q (pieces 22.6-30.2 Mb at CN 4) becomes a single
      7.5 Mb line; one line per event carried by several members (v1.7.0);
   c. v11: the 131 decided labels into training (karyotypes set apart: the
      X DUPs of an XXY received contradictory verdicts).
2. **CNN second opinion (1sexdecies)** — step 1 DONE (v1.3.0, 2026-11-04):
   `svengers tensors`, 13 channels x 256 points + 12 scalars per call
   (PASS/LowSQ >= 1 kb, ~1,200 per genome, ~20 s, ~1 MB). The cohort
   context is the panel of normals (median and MAD of 300 founders per
   bin), with no extra computation. Check: the PMS SHANK3 DEL (SQ
   0.38) reads as a clean loss (log2 -1.0, flat flanks, no heterozygous
   site, flat panel); a second rejected DEL from the same genome has the
   same profile. Remaining: cohort tensors, light 1D CNN (log1p on the
   counts), "RESCUE_CNN" use on the LowSQ calls once the IGV labels are
   available.
3. **BAF from the DeepVariant VCF, without a QUAL threshold** (finding of
   2026-11-04). The GHFC VAF track keeps `PASS && QUAL >= 30 && DP >= 10`;
   DeepVariant 1.9 gives QUAL < 30 to many true heterozygotes
   (mean VAF 0.47, 90% between 0.3 and 0.7). An ordinary genome loses
   ~40% of them; high-depth genomes (>= 60x) and all IP genomes
   up to 98% (het/hom ratio 0.05-0.85 instead of 1.33: 6 non-IP genomes
   out of 2,825, the 135 IP ones). Their BAF is almost empty: neither loss
   of heterozygosity to confirm a DEL, nor allelic imbalance for a
   DUP, and IBS0 relatedness reads parent-child pairs there as
   unrelated. Genomes at 11-18x lose sites through DP >= 10 but
   keep a normal het/hom ratio. The store reproduces the track
   identically: the import is not at fault. Plan: a command that rewrites
   `baf/` in the existing stores from the DeepVariant VCF (PASS, DP >= 10,
   no QUAL; no re-reading of the CRAMs, like `svengers wave`), bump of
   `SCHEMA_VERSION`. The model's BAF features change distribution:
   retrain (to do with v11 so as to retrain only once) and the usual
   deployment checks. Census: untracked script
   `runs/_sexpon/vaf_census.py`.
   Command `svengers baf` DONE (v1.8.0, store 0.6, old BAF kept as
   `baf_track/`, ~1 min per genome). On copies of 5 stores: autosomal
   het/hom 0.04 -> 0.91 (IP 75x), 0.08 -> 1.18 (91x), 1.40 -> 1.70
   (ordinary); karyotypes unchanged, x_het of IP males 0.08-0.12 ->
   0.02-0.03 ("single X" threshold 0.15). Remaining: run on the production
   stores, cohort check (karyotypes, IBS0 relatedness of IP parent-child
   pairs), then recomputation of the calls' BAFN/BAFDEV from the
   store (`svengers regenotype`, DONE in v1.10.0: reproduces the call's
   values on a non-rewritten store) and v11.
   Production stores rewritten on 2026-11-05 (2,973, 0 failures). Summary
   (`runs/_baf/`): relatedness — median KING of parent-child pairs 0.116 ->
   0.233 (IP -1.99 -> 0.205), IBS0 0.111 -> 0.018; 1,574/1,647 PO pairs of
   dataset_paper at KING >= 0.177 (none before), the other 73 mostly IP
   (BAF still incomplete: RefCall `./.` sites) or borderline; two pairs
   declared unrelated at KING 0.48 (duplicate or twins). Karyotypes:
   0 changes out of 2,989. v5n on re-genotyped calls: PASS per genome
   -1.1% DEL, -2.5% DUP outside IP (max ~3%). Mosaics: 9 low-fraction
   segments on wavy genomes passed the BAF test, hence the fraction
   floor of 0.3 on wavy genomes (DONE, v1.11.0); thereafter
   identical to the old BAF plus a 22q13 loss at f = 0.18 (to review in
   IGV). Separately: IP genomes keep ~22k non-genotyped RefCall sites
   (`./.`) at VAF 0.3-0.7 on chr20 (~6k on an ordinary genome),
   true heterozygotes rejected or paralogues: not taken up.
   Cohort genotyping on the new BAF: median de novo per child
   27 -> 42 on dataset_paper (v5n, run `runs/_v11/dp_v5n_rebaf`). Cause:
   at paralogous loci the new BAF has ~20 "het" per kb, each one a
   mismap penalty under a deletion; a common deletion genotyped
   0/1 in 98% of genomes drops to 40%, leaves CohortArtifact /
   HighCarrierFrac and produces de novo calls (16,800 rows over 1,503
   sites). chr12: 3,025 de novo (production), 3,653 (new BAF), 2,840
   (current code, original BAF). Genotyping therefore reads `baf_track` by
   default (`[genotype] rewritten_baf`). Paralogue-robust BAF likelihood
   (DONE, `[genotype] baf_paralog_p`): the hets of an event are ignored
   when a diploid stretch of that length would give as many with a
   probability < 0.05 (Poisson at twice the contig's het rate). chr12,
   new BAF: 3,653 -> 2,921 de novo (original BAF 2,840), violations
   280,264 -> 271,195 (268,078); no effect on the original BAF (2,862,
   267,883). Also counting 1 kb of flanks adds nothing (3,067 at 1e-2).
   The new BAF remains slightly worse than the old one for genotyping
   (diffuse low-QUAL hets, ~0.8% of positions on male chrX versus
   0.07%): stay on `baf_track` for genotyping, the new BAF serving
   mosaics, the karyotype and the model's features.
4. **v11**: IGV labels + transmission + GIAB, on the new BAF (point 3),
   usual deployment checks (PASS per genome vs production,
   recurrence, trios, HG002 out of sample). Train in it the integer CN read
   with the flanks on the edge bins (local branch `cn-partial`, events
   ≥ 1 kb): correct on chr22:49.32 Mb (0/1 -> CN1, nested -> CN0), but
   with v5n it loses 165 transmitted 1-5 kb PASS DELs over 20 trios for
   86 gained (PASS per genome -0.1%).
5. **GHFC export by the pipeline**: the Nextflow wrapper publishes directly
   to `samples|families|cohorts/…/svs/vedette/`, no more manual copying.
6. **Very noisy genomes at calling** (hidden-split budget, early
   detection, memory); complex **mosaics** and chrX; publication
   last.
7. **`CohortArtifact` unstable on a small cohort** (finding of
   2026-11-05). The flag fires on `AF > 0.02 and HWE_P < 1e-6`, exact HWE
   on the founders. On PMS (214 founders), ~140 common sites have a
   `HWE_P` stuck at the threshold (heterozygote excess typical of
   paralogues): a single founder fewer (a parent-child link added to the
   pedigree) doubles their p-value and 163 sites change FILTER (25 become
   PASS, AF 0.13-0.29). No effect on rare sites (identical family pages),
   but the status of common sites depends on the pedigree to within one
   individual.
   Leads: one-sided test for heterozygote excess (which is what signals a
   paralogue artifact; a heterozygote deficit is a different signal),
   threshold as a function of the number of founders, or heterozygosity F
   with a margin (|F| rather than p); check on dataset_paper (2,670) that
   the known artifacts remain flagged. The current status is in `HWE_P`
   (INFO, rounded to 6 decimals: round it less for auditing).
   Analysis (2026-11-05, `runs/_cohortart/`, untracked): the instability is
   not specific to the corrected link. Removing any founder flips
   ~150 sites on PMS (206 founders, 0.46% of the 31,928 flagged), 59 on
   dataset_paper (1,688 founders, 0.04%): it is the lot of any hard
   threshold, amplified by a small cohort, and the borderline sites are
   common (median AF 0.45), never on the family pages. A one-sided test
   for heterozygote excess or a threshold on F only reduces the flips by
   3-50% and would release the sites with a heterozygote deficit (557 PMS,
   4,420 dataset_paper: F 0.4-0.9, 88-93% < 1 kb, genotypes biased towards
   1/1): current rule kept. What a small cohort really lacks is
   power: of the 54,492 common PMS sites matched to dataset_paper, the
   latter flags 10,767 more (PMS alone: 76 more); 93% of the borderline
   PMS sites (p 1e-8..1e-4) are matched there, CohortArtifact at 97%.
   Lead: with `--known-sites`, take over the CohortArtifact of the
   reference site (power and stability of the large cohort). DONE (v1.10.0,
   `REF_ARTIFACT`, common sites only); PMS test
   (`runs/_test_refart_PMS`): 11,597 sites taken over, PASS 28,093 ->
   23,129 (all common), page rows -39%, PASS de novo candidates 2,203 ->
   1,302, no 22q13 deletion affected.

### State as of 2026-11-02 and priorities

In production: `models/v5n`, SQ >= 5, median PoN, sex chromosomes
corrected (point 6 below, DONE): `runs/dataset_paper` and `runs/PMS`
re-called, old versions archived as `runs/*_prod20261017`. dataset_paper
autosomes identical to the previous production; male chrX at the level of
the autosomes in transmission; chrY limited to the sites that the copy
number supports (`ChrYNoDepth`); de novo candidates of a sex-chromosome
aneuploidy kept out of the safeguard and out of the tables (an XXY boy
passed for a technical outlier, and a piece of his supernumerary X over
MECP2 for a de novo DUP).

Priorities:
1. **Swap follow-ups** (small): "events" pages that still show the
   pieces of the supernumerary X as de novo; 2 chrY PASS DUPs discordant
   father/son at 75%; a boy at chrX 0.84 / chrY 0.65 to review in IGV;
   old male stores in the archives to clean up after validation.
2. **Independent labels for depth**: targeted IGV batch of large
   depth DELs/DUPs, then CNN second opinion (1sexdecies), then v11 —
   no retraining on transmission labels alone (v9/v10 conclusion, still
   valid on the autosomes).
3. **Very noisy genomes at calling**: remainder of point 4 of the
   previous list (hidden-split budget, early detection, memory).
4. **Mosaics**: complex cases (point 3 below) and chrX.
5. Publication (plan below), not the priority.

### State as of 2026-10-14 (v1.1.0 released) and priorities

In production: `models/v5n` (v5 with `wave_score` neutralized, 1septdecies),
SQ >= 5, report v1.1.0. Two IGV reviews by Alexandre (138 + 44 cases,
`runs/validation/`, untracked) set the report rules:

- "below the NDD threshold" tier: true 1/10 in IP genomes (old, of lower
  quality), 1/13 in wavy genomes, 15/28 elsewhere → excluded for IP
  (`--low-quality-prefix`) and wave >= 0.128; at most 3 per genome;
- high-impact de novo: 56/63 true; the 7 false ones span 177 bp to 3.1 kb
  (including a 3.1 kb site in 3 families) → "< 1 kb" and
  "recurrent" alerts below 10 kb; de novo tables per event, no longer per
  site;
- breakpoints: fitted on the depth of the carriers against the
  non-carriers of the family (`refine_bounds`) — median start error of
  0.5 Mb brought down to a few kb on the 1st review; 2nd review: correct
  7/7 below wave 0.128, 9/18 above 0.2, the 2 false DUPs → deletions only,
  not above wave 0.2;
- mosaic fragments grouped under their segment; a deeper event inside
  stays separate.
- `svengers wave` filled in the `wave_score` of the old stores: 476/2,670
  wavy genomes on dataset_paper (97 visible before), ~100/319 on PMS.

Assisted review: HTML page per series with TSV export, carrier/family
depth plots (same 44 cases: 30 true in IGV, 38 with the plots — the plots
show the signal that the algorithm follows, they are not an independent
check in wavy genomes).

Priorities:
1. **Coverage waves at the source** — DONE: median PoN without SVD in
   production on PMS (2026-10-15) and dataset_paper (2026-10-16: sites
   361,359 -> 264,386, median de novo per child 194 -> 112, Mendelian
   violations 4.02% -> 3.15%, true large events from the IGV review PASS
   25 -> 28, false PASS 6 -> 4, CNDiscordant 58,588 -> 12). History: wavy
   genomes carry most of the remaining false positives in every review
   and made v6 to v8 fail. Panel of normals of unrelated founders
   (dataset_paper) + SVD denoising (already coded, never validated, cf. 4),
   tested first outside calling (wave_score before/after, 22q13 deletions
   intact), then PMS re-called with `--pon` and compared with v5n (PASS
   per genome by wave, sites, de novo, MVR). Sequencing metadata (batch,
   kit, date) requested for a per-batch normalization.
2. v9 and CNN second opinion (1sexdecies). **v9 not adopted (2026-10-17)**:
   v8 recipe on the calls with PoN (GIAB HG002/3/4 re-called with the PoN,
   24 wavy children labeled by transmission on dataset_paper PoN,
   wave_score neutralized). Wavy holdout much better (DEL P/R 0.945/0.647
   -> 0.957/0.810; DUP 0.262/0.049 -> 0.846/0.320), IGV review +8 true
   with no loss, DUPs >= 1 kb finally passing (31 -> 976 over 40 children,
   69% transmitted); but dataset_paper refiltered in full: PASS de novo
   (GQ >= 20) median per child 27 -> 63 (1-10 kb 14 -> 44, >= 10 kb
   2 -> 5), MVR 3.15 -> 3.22%, transmission 0.894 -> 0.878, and HG002 out
   of sample F1 DEL 0.451 -> 0.405. v9's labels do not see the failure:
   the apparent 1-10 kb depth de novo calls, almost absent from the labels
   (the GIABs have no clean parents apart from HG002, the wavy children
   are few). v10 lead: cohort labels on hundreds of dataset_paper PoN
   trios (de novo with confident parents = FP, transmitted = TP,
   excluding the evaluation trios), evaluation by the cohort trio metrics
   (now ~1 h for 2,670 genomes).
   **v10 not adopted (2026-10-17)**: cohort labels on 208
   dataset_paper PoN trios (non-IP parents, wave < 0.128; TP = carrier
   parent WITH its own call, `label-trio --parent-calls`; FP = two
   confident 0/0 parents without a call; common cohort sites and calls
   present in >= 10 of 100 unrelated founders discarded: 1.68 M -> 141 k
   labels) + GIAB PoN. Panel of 40 trios: median apparent de novo per child
   526 (v5n) -> 912, PASS DUPs >= 100 kb 0 -> 13 per child, noisy genome
   1,953 -> 4,526 PASS; HG002 out of sample F1 DEL 0.451 -> 0.412 (CMRG
   0.452 -> 0.350). Cause: the transmission labels are biased for depth
   calls — an FP requires two confident 0/0 parents, rare in noisy
   regions, a TP a single parent with a call: 75% of the large depth DUPs
   labeled TP; v10 passes 712 depth DUPs >= 10 kb over 40 children
   (v5n 7), of which 291 without any parental call. **v9/v10
   conclusion**: without independent truth for depth calls, retraining
   shifts the errors; v5n + PoN remains the best measured compromise.
   Next for the model: independent labels for depth (targeted IGV batch of
   large depth DELs/DUPs, CNN second opinion), no new retraining on
   transmission labels alone.
3. Complex mosaics: mosaic cut by a deeper deletion
   (chr13), gain/loss rearrangement in 7q11.23.
4. **Very noisy genomes at calling** — partly DONE (2026-10-19): the
   search for the SR pile of each PE cluster went through all the piles of
   the contig (quadratic; 1.9 M piles on chr20 of a noisy genome) ->
   binary search, identical outputs; whole noisy genome 8 h 56 -> 3 h 10,
   identical VCF; memory peak ~48 GB -> Nextflow retry at 72 GB. Remaining:
   the hidden-split budget exhausted on every contig (~2 h of the
   remaining time), and the memory peak after the scan. Initial finding:
   (noted during the dataset_paper re-call with PoN, 2026-10-16): 5 genomes
   from the same batch killed by OOM at 24 GB in ~5 min, then 5 to 9 h per
   call instead of ~1 h (~5 M weak pileups per chromosome, hidden-split
   alignment budget exhausted on every contig). They delay the whole
   cohort by ~8 h. To do: detect the profile early (pileup density, raw
   wave) and cut the hidden-split search earlier or bound it in time;
   reserve more memory for these genomes from the start rather than
   waiting for the failure; flag them in the report as genomes to treat
   with caution.
5. **Sex-chromosome karyotype** — DONE (2026-10-17, `karyotype.py`, "Sex
   chromosomes" section of the report; 2,989 GHFC genomes: 1 XXY, 1 XYY,
   1 XXX, 1 probable X mosaic, 8 partial Y losses, 4 males with a
   heterozygous X = contamination?, 5 PMS families with the parents' sexes
   swapped). Initial plan: the mosaic scan covers
   only the autosomes, a whole-chrX gain (probable XXX on dataset_paper)
   is seen only through the calls. Explicit X/Y check per genome: mean CN
   of chrX and chrY (and BAF of chrX outside PAR) against the pedigree
   sex → XXY, X0, XXX, XYY and sex-chromosome mosaics, sex
   discrepancies; a row in the table of constitutional chromosomal
   anomalies. Small job, doable in parallel with 2.
6. **Male chrX/chrY with the PoN** — DONE (2026-11-02): males
   re-called (1,647), cohorts redone and in production. Along the way, the
   depth noise (sigma0) was fitted in a 0.7-1.3 window: on the chrX of a
   male at 0.5, ~170 bins out of 147 k, sigma0 0.02 instead of 0.08, and
   every one-bin dip became a hemizygous deletion — fixed
   (window around ploidy/2). Sex-aware cohort steps (HWE on the
   females, carrier fraction per sex, X transmission mother → son, Y
   father → son). chrY: 453 of the 525 PASS sites were small junction-only
   calls in the Y repeats, carried by up to two thirds of the males at
   their baseline CN, father/son discordant at 40-54% ->
   `ChrYNoDepth`; 72 sites remain (DELs discordant at 20%, chr7 36%).
   History: `PanelOfNormals.normalize` divided the chrX
   of a male (~0.5) by the median of the panel's males (~0.5): ratio 1.0
   against a baseline at 0.5, the whole chrX (and chrY) called DUP CN 2
   in the 1,501 males of dataset_paper (median 156 Mb) and those of PMS
   since the switch to the PoN. v5n rejects them (SQ ~0.01), but cohort
   genotyping reads every male as a carrier of the chrX DUP sites (8,637 of
   the 10,182 filtered chrX sites have > 50% male carriers,
   CohortArtifact): a real chrX DUP in a boy (MECP2) risks being discarded.
   Fix: ratio brought back to the ploidy/2 level; male genome re-called:
   identical autosomes, chrX depth DUPs 156 -> 0.2 Mb; female genome
   identical. Re-call the males (dataset_paper 1,503, PMS 153)
   then redo merge -> report.
   **IGV review v4 (2026-10-30, 55 depth calls >= 10 kb, 40
   children, blinded v5n/v10)**: DEL v10 only 16/16 true, v5n only 1
   true / 3 false (DUSP22) / 2 uncertain, both 8/8; DUP both 4/4
   (LARGE_DUP_BAF rescue), v10 only autosomes 3 true + 1 uncertain,
   v10 only chrX/chrY 0/17 — all of them this bug (normal XY males). v10's
   transmission labels were contaminated by it (every father a "carrier"
   of his whole chrX). To redo after the re-call: cohort labels then v11.
7. dataset_paper paper: out-of-sample figures, comparison with the
   other callers on the trio metrics, figures.
### Publication plan (2026-10-14) — not the priority

Alexandre's decision: first a top-notch tool, publication afterwards. What
submission will require, from blockers to finishing touches:

1. **Out-of-sample evaluation (blocker).** The headline figure (HG002, P 0.95)
   is in-sample: HG002 is in v5's training set. Main figure with the
   production model retrained without HG002; additional public truth sets:
   HGSVC samples (assembly-derived variants) with their public Illumina 30x
   data (1000 Genomes, NYGC), HG002 Q100.
2. **Comparisons (blocker).** Depth-based CNV tools in addition to SV
   callers: CNVnator, GATK-gCNV, Canvas (and DRAGEN if accessible) alongside
   Manta, DELLY, dysgu; same data, default settings; on GIAB (Truvari) and
   on a cohort (de novo rate, Mendelian violations).
3. **Orthogonal validation.** The ~180 IGV validations (runs/validation,
   untracked); results from another technique (array, ACPA, MLPA, clinical
   reports), in particular the known SHANK3 deletions in PMS, for a
   real-world sensitivity.
4. **Coverage waves.** A method (panel of normals + SVD, in progress)
   evaluated before/after, rather than the current neutralization and
   exclusions.
5. **Framing.** DUPs clearly weaker than DELs: present them as they are or
   centre the paper on deletions. Ablations (evidence types, learned filter,
   cohort re-genotyping). Mosaics: simulated fractions and confirmed cases.
   Costs (CPU and storage per genome, 2,670 genomes).
6. **Reproducibility outside GHFC (editor blocker).** Generic config and a
   downloadable GRCh38 resource bundle (masks, GC, popCovStdev); model
   trained on public data only (v5 uses 16 private GHFC children) or a data
   access statement; public benchmark pipeline that regenerates the figures;
   archived release with a DOI (Zenodo), public container; Bioconda/PyPI
   only with Alexandre's agreement.
7. **Writing.** Application note (fewer comparisons) or full paper: to be
   decided; methods, figures, cohort ethics.

### 0. (archive) Production: models/v5 + threshold SQ>=5 (2026-09-18)
v5 = GIAB + 16 GHFC wavy-genome trios (transmission labels, `svengers label-trio`).
Holdout of 8 wavy-genome children: DUP PR-AUC 0.908 (v4) -> 0.988; DEL 0.989
-> 0.997.
HG002: no regression (SQ5: P 0.951 / R 0.334 vs v4 0.953 / 0.337).
Optimized filter DONE: batched scoring (2026-09-18, 6h05 -> 21 min) then
feature extraction with a per-contig depth cache (2026-09-21, 21 min ->
13 s/genome with v5; VCF and feature TSVs identical to the reference).

### 0bis. (archive) Production v4 (2026-09-17)
GHFC pipelines swapped to models/v4 (30 features, wave_score),
low_sq_threshold=5 in configs/ghfc.toml (shifted calibration: SQ5 v4 =
SQ7 v3 operating point; HG002 v5.0q: P 0.953 / R 0.337 / F1 0.498).
Next expected gain: retraining enriched in GHFC wavy genomes.

### 1. ML filter training — v1 DONE (2026-09-09)
`models/v1` models trained on HG002 40x NovaSeq PCR-free aligned GHFC-style,
GIAB v5.0q stvar truth, Truvari labels (5644 TP / 60373 FP).
- DEL: PR-AUC(OOF) 0.952. Genome-wide eval vs truth at SQ≥7: **precision
  99.1% (FDR 0.9%, gnomAD 2.8% target beaten ×3), recall 25.9% overall**;
  by size: 300 bp-1 kb 53%/99.6%, 1-10 kb 66%/99.2%, 10-100 kb 48%/100%.
  Weaknesses: 50-300 bp (13.7%) and segdup (6.2% recall at 98.3% precision).
  Raw candidate recall (pre-filter): 40.8% ≈ Manta.
- DUP: PR-AUC(OOF) 0.643, SQ≥7 → FDR 7.8% ≈ 7.3% target. Vs INS truth:
  precision 97.2%, recall 1.4% (the INS truth mostly contains non-DUP
  insertions — a clean DUP eval requires isolating the tandem-dup subset of
  the truth, see the bench backlog).
- Backlog v1.1: (a) recall of small DELs 50-300 bp (SR-only junctions pruned
  too hard, truvari matching of DELs in sequence representation); (b) DUP FPs
  upstream (stricter everted nomination, both sides required); (c) HG003/4
  semi-supervised; (d) tandem-only DUP eval DONE (12321/28942 tandem INS;
  SQ≥7: 91.3% precision vs Manta 60.1% at equal recall, raw reservoir 5×
  Manta); (e) rerun with a whole-genome founders PoN; (f) normalize the
  SVLEN/END convention of the VCFs (htslib forces END=POS+|SVLEN| on
  symbolic alleles, +1 vs our half-open — consistent everywhere but to be
  cleaned up before publication); (g) retrain v4 with the wave_score
  feature; (h) merge rescue criterion adaptive to cohort size (2670 samples:
  80% of sites are RESCUE with ΣSQ>=8/2 carriers — tighten to a carrier
  fraction + mean ΣSQ, or calibrate on the CohortArtifact artifact rate);
  (i) family workspaces in a single batched pass (the per-family VCF subset
  re-reads the whole cohort VCF: 1207 families x 635M genotypes = infeasible
  — assemble from the per-sample gt.tsv files, each file read only once);
  (j) report: assemble the genotype matrices from genotypes/*.gt.tsv rather
  than re-decoding the cohort VCF (1.4 billion genotypes ≈ 1 h pysam vs
  ~5 min TSV).

### 1undecies. Breakpoints of family events: the carrier's, not the representative's (2026-09-30)

Alexandre, on FAM-10 (SHANK3 deletion): IGV places the deletion at
chr22:45,145,667-50,808,468, the family page displayed 44,094,543-50,800,002.

The per-sample caller was right: the proband's own call starts at
**45,146,001**, 300 bp from IGV. Two errors came afterwards, in the report:

1. **Cluster representative coordinates.** At merge, the proband's call
   (45.15-48.99 Mb) joined a site whose representative — another PMS
   patient, the cohort gathers 22q13 deletions of all sizes — starts at
   44,116,001. Re-genotyping gives the proband 0/1 GQ 99 (79% of the site is
   indeed deleted in them), and the row inherited the other patient's
   breakpoints: 1 Mb too much. Same defect the same day on **FAM-09** (IGV
   43,784,666, page 43,570,001), but there the own call comes back in nine
   pieces, none of which overlaps the site at 50%: a "one own call" rule is
   not enough.

   Final rule (`own_span`): the carriers' own calls that overlap the site
   (≥ 10% of their length inside the site) are chained with the collapse
   rule; the best chain replaces the site's breakpoints if it overlaps it at
   RO ≥ 0.5. The chain's breakpoints come only from calls that are **PASS,
   not nested in a longer call, and whose cohort site is PASS**; the other
   own calls only serve as links. Each guard answers a real PMS case:
   - filtered calls as breakpoints → 57,656 rows moved on the strength of
     rejected calls (hence "PASS only");
   - FAM-08: the middle (138 kb) of the deletion is a LowSQ call; without
     linking through filtered calls, the event was cut in two;
   - FAM-11: the duplication joined a site starting 141 kb inside it; a
     "never beyond the site" rule cut it in two (so the chain may enlarge
     the row);
   - FAM-14: a 23 kb subtelomeric PASS call, whose site is carried by all
     319 samples, pushed the end beyond chr22 (hence "cohort site PASS"); a
     call touching the site over only 1 bp did the same (hence the 10%);
   - FAM-06: five PASS calls of 0.1-0.8 kb nested in a 712 kb LowSQ
     duplication became breakpoints and cut a 5.6 Mb duplication (hence
     "not nested").

   Result on the 120 PMS families, against the 30/09 12:24 report: FAM-09
   moves to 43,785,001 (335 bp from IGV), FAM-10 / FAM-11 / FAM-15 /
   FAM-01 / FAM-08 / FAM-14 / FAM-06 unchanged, no row beyond the end of a
   chromosome; 3% of rows move (median 1.6% of their size); among events
   ≥ 1 Mb, six change breakpoint, each towards the carrier's own PASS call
   (FAM-02, FAM-07, FAM-09, FAM-12, FAM-13, FAM-05 by 1 bp).
2. **Bridge over a small neighbor.** A 0.6 kb polymorphism (NCAR 38, mother
   at GQ 11 hence non-carrier) 21 kb from the edge was glued on by the gap
   tolerance (2% of 4.9 Mb = 97 kb). Fix: crossing a gap requires the
   smaller of the two neighbors to be at least as long as the gap; the
   fragments of a single event are large, neighboring polymorphisms are not.

End of the event: 50,800,002 versus 50,808,468 in IGV. This is not a
breakpoint discrepancy: the 50,800-50,819 kb segment is subtelomeric
(LowMap, MaxN), carried by all 319 samples, hence non-PASS. 50,808,468 is
the chromosome end minus the 10 kb telomeric gap.

Still open, this time in the genotyper: a partial carrier gets GQ 99 on a
site that extends 21% beyond their event. A CN test on the site edges
(first and last 10%) would flag this partial overlap for carriers
**without an own PASS call**, which the report fix cannot cover.

### 1septdecies. `wave_score`: v5 floods wavy genomes; v5n; v8 not retained (2026-10-08)

**Production tracks reproduced** for GIAB: `by1000.bedgraph` =
`mosdepth -n -m -x -Q 40 --by 1000` (median, fast mode, MAPQ ≥ 40:
64,445/64,445 identical bins on chr20 of a GHFC genome); `vaf.bedgraph`
= DeepVariant 1.9.0 (`run_deepvariant --model_type WGS`, parameters = those
of `ghfc-ngs`), PASS sites, QUAL ≥ 30, DP ≥ 10 (0 extra, 0 missing).

**Discovery along the way: v5 learned a step on `wave_score`** (per-genome
feature) between 0.10 and 0.13. Same calls, same features, only the wave
changes: PMS genome 4,794 → 29,444 calls at SQ ≥ 5 (wave 0.10 → 0.13);
HG002 production tracks (wave 0.164) 42,082; a dataset_paper genome forced
to 0.16: 5,070 → 42,764. Production affected:
- PMS (schema 0.4 stores, wave present): 102/319 genomes ≥ 0.12, median
  PASS 12,600 versus 5,500;
- dataset_paper: 2,517 stores in schema 0.3 (wave absent = −1, hence a
  "healthy" v5 by accident), but 153 in 0.4 of which **97 ≥ 0.12** — very
  probably the ~81 "load ≥ 5×" genomes of the reports.
HG002 with production tracks, v5: DEL P 0.456 / R 0.327.

**v5n** (`models/v5n`, commit e4c58b1) = v5 estimators unchanged,
`wave_score` neutralized to −1 at scoring (`neutral_features` in the
manifest; `svengers train --neutralize`). Identical to v5 on any genome
without a wave or below the step. PASS per genome: normal genomes
identical; wavy genomes 10-18k → 3.3-6k (junction DUPs < 1 kb 10-11k →
0.4-0.9k). HG002 production tracks: DEL P 0.859 / R 0.367 / F1 0.515 (v5:
0.381); CMRG DEL F1 0.563.

**v8 not retained.** v5-like + GIAB on production tracks + cleaned wavy
labels + `wave_score` neutralized. Compared with v5n:
- PASS per genome 2.8-5.4k (v5n 3.3-5.6k); junction DUPs < 1 kb divided
  by ~2; depth DELs ≥ 1 kb stable (not the v7 explosion), except for one
  wavy PMS genome (457 versus 119, 17% recurrent);
- cleaned wavy holdout: DEL P 0.970 / R 0.585 (v5n 0.829 / 0.635); DUP
  P 0.802 / R 0.163 (v5n 0.179 / 0.831);
- HG002 out of sample, production tracks: DEL P 0.898 / R 0.290 /
  F1 0.439 versus v5n 0.515; CMRG DEL F1 0.472 versus 0.563.
More precise, but the recall lost on out-of-sample HG002 and CMRG
outweighs it. To dig into before a v9: the share due to neutralizing the
wave (which perhaps carried a real signal for wavy genomes) and the share
due to the redone GIAB labels.

### 1sexdecies. CNN second opinion on the store tracks — PLAN (2026-10-06)

Alexandre's question: do for DEL/DUP what DeepVariant does for
SNVs/indels, i.e. classify an "image" of the event.

**Starting observation.** The filter (gradient boosting, `filter/schema.py`)
only sees summaries: `dhffc`, `depth_cv`, `frac_consistent`,
`n_bins_informative`, `flank_asym`… never the *shape* of the signal.
Textbook case (PMS, 2026-10-06): 9 kb deletion of the 3' end of SHANK3 in a
proband, CN 1, RD 0.50, GQ 99, unique in the cohort, confirmed in IGV in
one second (two clean steps at 0.5, flat flanks at 1.0) — SQ 0.38 for a
threshold of 5: depth only, no heterozygous site (GC-rich 22q13.33),
`depth_cv` 0.30. Conversely a GC coverage wave has an undulating shape
without a step that summaries distinguish poorly (cf. the ~70 CN ~1.83
waves of the dataset_paper mosaic scan, 1quaterdecies). The "below
threshold, NDD genes" tier of the family pages (2026-10-06) catches this
case for review; the CNN would aim for the score itself to see it.

**Literature.** Samplot-ML (Belyeu et al., GigaScience 2021: CNN on samplot
images of DELs, false-positive filtering); Cue (Popic et al., Nat Methods
2023: multi-channel images of alignment signals, object detection);
CNV-espresso (exome: depth images of the sample *and* of reference samples,
in silico confirmation); DeepSV; SVision (long reads).

**What we do NOT do: render PNGs.** The signal of a CNV of a few kb is not
base-resolved like that of an SNV: it is already condensed in the store
(`depth`, `baf`, `pe`, `sr`, `ins`, `cig`). The store *is* the image;
rasterizing loses information and is expensive. Invariant respected: the
second opinion reads only the store, never the CRAM.

**Plan.**
1. *Tensors* (`svengers/filter/tensor.py`, new): for each call,
   window = call + a flank of the same length on each side (bounded, min
   10 kb), resampled to L = 256 points; channels:
   - log2 of the diploid-normalized depth ratio, and mask;
   - BAF: density of heterozygous sites and median |BAF − 0.5| per point;
   - discordant pairs (oriented DEL/DUP), split reads, insert-size shift
     (`ins`), CIGAR indels;
   - **cohort context**: median and MAD of the depth ratio of the cohort's
     genomes (or of a reference panel) over the same bins — the lever
     against systematic artifacts (chr16:28-31 Mb dip across the whole PMS
     cohort, GC coverage waves). Cost: one per-bin summary computed once per
     cohort (like `svengers mosaic`, ~1-4 min), not per call;
   - scalar channels concatenated after the convolution: log10 size, type,
     source, `popcov_stdev`, `segdup_ovl`, `lcr_ovl`.
   Extraction ~ms per call; to be restricted to pre-filtered calls (hard
   filters passed), not to the ~35k raw calls per genome.
2. *Model*: lightweight 1D CNN (3-4 conv + pooling blocks, < 1M parameters),
   one per type or a shared one with the type as input; CPU sufficient for
   inference. Output = calibrated `SQ_CNN` (isotonic), written in FORMAT
   next to `SQ`. The GBM feature schema is not touched (append-only).
3. *Labels* — the hard part, cf. v6/v7 (1quindecies):
   - GIAB HG002/3/4 truth (Truvari), BUT with depth produced as in
     production (GHFC-style 1 kb bedgraph) — prerequisite shared with v8;
     otherwise we replay v7's CRAM/bedgraph mismatch;
   - cleaned trio transmission (`label-trio --known-sites`);
   - **IGV validations from the family pages** (true/false/uncertain, TSV
     export): human labels on exactly the hard cases. Goal before training:
     a few hundred, including "below threshold" calls.
   - augmentation: left/right symmetry, slight breakpoint jitter, depth
     noise; no injected synthetic events as the sole source.
4. *Two-stage use*:
   a. **second opinion on LowSQ calls only** (like SHANK3): a call rejected
      by SQ but with a high SQ_CNN becomes "RESCUE_CNN", with a limited
      scope;
   b. if (a) holds: SQ/SQ_CNN combination (stacking) for all calls.
5. *Deployment criteria* (v6 lesson, "pre-deployment check" memory): PASS
   per genome **by type, source and size** vs v5 on ≥ 5 genomes including
   wavy genomes; recurrence of the new PASS calls in unrelated individuals;
   HG002 out of sample; cleaned holdout; and, for (a), number of rescues per
   genome (target ≲ 1) checked by IGV review on a sample.

**Estimate**: tensors + cohort context ~2 d; training/evaluation ~2 d once
the IGV labels are available and v8 (production-style depth for GIAB) is
done. Do not start before these two prerequisites.

### 1quindecies. Training/production mismatch on `gq` — to be fixed by v6 (2026-10-05)

HG002 benchmark of v1.0.0 (docs/benchmark_hg002_v1.tsv) against that of the
v5 model (09-18): DEL 50-300 bp recall 0.225 → 0.210, LCR 0.192 → 0.172,
overall F1 0.495 → 0.483; CMRG identical. Autopsy of the 132 lost true DELs
of 50-300 bp: 103 are still called (JUNCTION) but fall into LowSQ, 21 are
missing. The filter is not to blame (old calls refiltered with the v1 code:
same median SQ 5.16). Between the two call sets, **a single feature changes
for the 104 matched candidates: `gq` (+10 at the median)**; putting the old
`gq` back on the v1 features makes 95 of the 104 PASS (median SQ 4.35 →
5.16).

Origin: `384b25a` (09-09 17:34) replaces the Poisson likelihood of the
junction read count with a negative binomial. The GIAB HG002/3/4 calls used
to train v3 → v5 date from the morning of 09-09 (Poisson); all production
calls (GHFC cohorts, HG002 v1) are later (negative binomial). The v5 model
is therefore serving in production a feature whose distribution has changed
since its training. The v5 benchmark (on the 09-09 calls) slightly
overestimated production; the v1 one is the faithful one.

**v6 trained on 2026-10-05** (`models/v6`, not yet in production): same
data and labels as v5, GIAB calls redone with the v1 code, features with a
live `lcr_ovl`. Holdout of 8 wavy-genome children (production features,
outside the training of both models), at threshold SQ ≥ 5: DEL P 0.974 /
R 0.714 (v5) → P 0.978 / R 0.973 (v6); DUP P 0.795 / R 0.728 → P 0.878 /
R 0.960 — the `gq` mismatch was therefore costing ~25 recall points in
production. HG002 with a v6 model trained WITHOUT HG002 (first
out-of-sample figure: v3-v5 were evaluated on an HG002 seen in training):
DEL R 0.365 / P 0.870 / F1 0.514 (v5 in-sample: 0.323 / 0.955 / 0.483;
Manta 0.376 / 0.934 / 0.536; dysgu 0.389 / 0.879 / 0.539); best F1 of the
table from 300 bp up (300 bp-1 kb 0.722, 1-10 kb 0.831, 10-100 kb 0.706);
CMRG F1 0.573 (v5 0.531). docs/benchmark_hg002_v6.tsv.

**v6 withdrawn from production on 2026-10-06** (put into production on
05/10, rolled back before any report). At threshold SQ ≥ 5, v6 lets through
12-17k calls per dataset_paper genome versus 4.7-5.6k with v5 — mostly
small junction DUPs < 1 kb (5-8.6k per genome) — and the merge goes from
523k to 1.92M sites (cohort-filter killed at 16 GB). The wavy holdout
(P 0.878 / R 0.960 for DUP) was misleading: it only scores calls with a
clean transmission label and discards the ~25k "ambiguous" calls per child,
which is where these DUPs are. Out-of-sample HG002 already said so (246 DUP
FPs versus 46, DEL precision 0.870).
Lesson: a model is validated on the NUMBER of PASS calls per genome and on
the cohort's trio metrics (MVR, apparent de novo) before production, not
only on a holdout with filtered labels.
To do for a usable v6: (1) clean the transmission labels (a recurrent
artifact present in a parent is "transmitted"; exclude common sites via
`--known-sites`, cf. 1ter); (2) calibrate the threshold on a target
precision measured on the full set of calls (out-of-sample HG002 + apparent
de novo rate), not on the clean labels; (3) deployment safeguard: compare
the number of PASS calls per genome with the model in place.

**v7 (2026-10-06), not retained either.** v6 + cleaned transmission labels
(`label-trio --known-sites`: no label on common sites or sites flagged as
artifact; 99% of the small "TP" junction DUPs fell there).
Assessment on four criteria:
- PASS per genome (6 dataset_paper genomes): normal genomes +12 to +26%
  vs v5 (v6: ×2.5 to ×3.4); overloaded wavy genome −40%; other wavy genome
  +66%.
- wavy holdout (cleaned labels, rare/moderate sites): DEL P 0.938 / R 0.870
  (v5 0.829 / 0.635); rare DUPs weak for all (v5 P 0.18, v7 0.27 with
  recall halved).
- HG002 out of sample: DEL P 0.883 / R 0.358 / F1 0.509 (v6 0.870 /
  0.365 / 0.514); in-sample v7 0.945 / 0.345 / 0.506 vs v5 0.955 /
  0.323 / 0.483; DUP (INS space) precision 0.939 vs 0.885; CMRG F1 0.573.
- BUT the extra PASS calls on cohort genomes are **depth-only DELs ≥ 1 kb
  ×6 to ×18** (269 → 1679; 169 → 3142 in a wavy genome), with no junction
  support, 45% present in the 3 unrelated individuals tested (v5: 28%) —
  very probably recurrent depth artifacts. Invisible on HG002: GIAB depth
  comes from the CRAM (giab.toml), cohort depth from the precomputed 1 kb
  bedgraph. Once the wavy labels (the only "bedgraph" examples) are cleaned,
  the model learns depth on CRAMs and serves it on bedgraphs: a new
  training/production mismatch, on the depth features this time.
- The transmission test is not conclusive (90% present in a parent), for
  the same reason as the labels: a recurrent artifact is "inherited".

For a usable model (v8): produce for HG002/3/4 the same depth as
production (GHFC-style 1 kb bedgraph, call on the fast path), redo features
and labels, retrain; deployment criteria = PASS per genome **by type and by
source** vs v5, recurrence in unrelated individuals, HG002 out of sample,
cleaned holdout.

Fix (done): **v6** — re-call HG003/HG004 with the v1 code (HG002 already
is, `resources/giab/hg002_v1/`), recompute the features, same labels, and
retrain (with `lcr_ovl` now live and realistic large events, cf. 1ter /
1nonies); re-benchmark. Safeguard to add: the model manifest should store
reference quantiles per feature, and `svengers filter` should warn when the
served distribution departs from them (drift test) — this mismatch would
have surfaced at the first production run.
### 1quaterdecies. Chromosomal mosaics (`svengers mosaic`, 2026-10-02)

Checked by Alexandre in IGV on father SMP-02 (family FAM-03): a single
deletion chr10:83.36 Mb-end where SVengers reported 57 CN=1 fragments. Depth
gives a **homogeneous CN ≈ 1.33** over 50 Mb: a loss in ~2/3 of the cells. The
same clone (fraction 0.62-0.67) carries 12 segments ≥ 5 Mb — losses of 2q, 6q,
chr9, 10q, 12p, 17p, gains of 4p and 17q (isochromosome 17q profile) — with
the expected allelic imbalance (|BAF − 0.5| 0.20 on 10q, 0.05 in the mother;
expected 0.25 for a loss, 0.125 for a gain at f = 0.67). Not transmitted to
the children: a somatic clone (culture artifact of a lymphoblastoid cell line,
or a clonal process in the sample), not a constitutional event — flagged to
Alexandre as a possible incidental finding.

Why SVengers missed it: the CN-HMM only has integer states; a ratio of 0.665
alternates between CN 1 and CN 2 (fragments), and gains at CN 2.67 fall on no
state (not called).

Pass added (`src/svengers/mosaic.py`, `svengers mosaic --run-dir`, ~1 s per
sample from the store): 1 Mb windows, binary segmentation, cell fraction
under the single-copy-change hypothesis (loss CN = 2 − f, gain CN = 2 + f),
then two filters: **homogeneity** (≥ 75% of windows within ±0.2 of the
segment CN — a mixture of integer levels is not a clone) and **BAF** (excess
|BAF − 0.5| ≥ 40% of the expected value for f). Mosaic if 0.15 ≤ f ≤ 0.85;
segments at f ≈ 1 (large constitutional gains, e.g. the 8q24 duplications of
the PMS translocations) are written but without an alert. The report reads
`samples/<S>/<S>.mosaic.tsv`: alert "mosaic 67% chr10:83-133 Mb" on rows that
fall within a segment of the carrier, "mosaic carrier" on the carrier's other
rows, "Chromosomal mosaics" table.

Calibration of the homogeneity filter on PMS: without it, a constitutional
3 Mb 16p13 duplication next to a deleted window (per-window CN 2.9 / 2.5 /
3.0 / 1.15) and a 22q13 deletion at an intermediate edge (1.66 / 1.65 / 1.03 /
1.02) passed as mosaics; with it, only one proband remains, at a homogeneous
CN 1.60-1.66 over 22q13 (43-50 Mb) with the BAF expected for f = 0.35 — a
candidate mosaic 22q13 deletion (ring 22?), to be checked.

Threshold lowered to 3 windows (3 Mb) the same day: Alexandre saw in IGV a
mosaic terminal 22q13 deletion (47.7 Mb-end, CN 1.62, BAF consistent with
f = 0.38) that the 5 Mb threshold missed, as did the main caller (no call: CN
too close to 2 for the HMM). The last partial window of the contig is kept
and the end bounded to its length. Effect: PMS 1 → 3 carriers, all terminal
22q13 mosaics (f 0.35-0.56); dataset_paper 3 → 4, the new one being a 3 Mb
7q11.23 DUP at CN 2.68 (segdup-rich region, mosaic or underestimated
constitutional duplication — to be looked at in IGV).

Limits: autosomes only (X/Y ploidy absent from 0.3 stores); events < 3 Mb out
of reach; complete deletions (f ≈ 1) no longer have heterozygous sites and
are not confirmed by the BAF (that is the main caller's job). To do:
integrate the pass into the Nextflow pipeline, and let the caller benefit from
the fraction (fractional states, or merging the fragments of a mosaic segment
into a single call).

### 1terdecies. Carrier QC flag on family rows (2026-10-01)

Request from Alexandre after the dataset_paper report: the impact score only
looks at gene content, and family FAM-03 took up about ten of the 57 rows
≥ +1.00 (PTEN, FGFR2…) with children at 478 and 94 de novo candidates.
"Carrier QC" column + filter (all / no alert / with alert): wave, PASS site
burden ≥ 5× the median, de novo safeguard.

Two defects found while validating it on PMS (49 of the 58 rows ≥ +1.00
initially flagged, whereas these are the true SHANK3 deletions):

1. **The absolute-threshold de novo safeguard flagged everyone**: > 20
   candidates = 664 of the 665 dataset_paper children and the 77 PMS
   children, and their de novo tables were hidden in the main report.
   Threshold now relative, like the carrier burden: 5× the cohort median
   (dataset_paper median 75, p90 367, p99 3173 → 63 children, including
   SMP-03 of family FAM-03; PMS median 100 → 1 child), the floor of 20
   remaining the threshold for family runs (< 20 children). Counting is done
   in events (sites merged by the collapse rule); my hypothesis that a large
   de novo counted for dozens of sites was wrong, the effect is marginal
   (SMP-01: 116 sites, 114 events).
2. **The wave score is inflated by large CNVs**: 54% of carriers of a chr22
   DEL ≥ 1 Mb are "wavy" against 13% of parents; chr22 alone is at 0.4-0.6 in
   them (low-frequency energy of a coverage step). Excluding the contigs that
   carry a large own call does not work (all 319 PMS samples have one, often
   artifactual, and the signal of the truly wavy genomes collapses). Chosen:
   on each row, the carrier's wave is the mean of their per-contig scores
   **without the row's contig** (`row_wave_score`) — 26 → 5 carriers flagged
   on their chr22 DEL; a truly wavy genome oscillates everywhere and keeps its
   flag.

4th alert "large CNVs" (2026-10-02): at least 3 events ≥ 1 Mb with a PASS own
call in the carrier (`big_cnv_counts`, per family page). On dataset_paper,
family FAM-03 was flagged by none of the three alerts: its 13 rows ≥ 1 Mb are
carried by father SMP-02 (wave 0.112 recomputed from his 0.3 store, normal
burden, no child), 13 non-transmitted deletions on 5 chromosomes — cohort
median 0, p99 1, max 13; the threshold of 3 flags 2 samples out of 2673
(SMP-02, SMP-04). Complex rearrangement or sample problem: to be settled by
IGV.

PMS result: rows ≥ +1.00 with a PASS own call, 51 without alert / 7 with
(FAM-09 and FAM-06 carry two large CNVs, chr22 DEL + chr8 DUP from a
translocation; FAM-14 wave 0.184). Still to do: the wave score itself (stored
at call time) should exclude called CNV segments — to be recalibrated against
the wave-detector labels (docs/wave_qc_validation.md).

### 1duodecies. Two points raised in the dataset_paper report of 2026-09-30

- **Call ends beyond the chromosome (caller, pre-existing) — fixed in the
  code, data to be redone.** 24 sites of `cohort/sites.vcf.gz` end up to 1 kb
  past the end of their contig (chrX 156,041,002 for 156,040,895; chr12
  133,276,002 for 133,275,309): the last 1 kb bin of a depth segment was not
  bounded to the contig length. 1362 rows of the family pages inherit it
  (identical before and after `78fe242`). Fixes: `f8174eb` bounds the end at
  write time (`vcfout.clamp_end`, shared by `write_vcf` and `run_merge`;
  internal bound `length - 1` because htslib writes END = POS + |SVLEN|, cf.
  (f): the written END falls on the last base of the contig); `a0d65be`
  bounds the own calls read by the report (`own_call_index`);
  `_read_cohort_vcf` also bounds the sites, for rows without an own call.
  Remaining: re-merge of dataset_paper and PMS (the published `sites.vcf.gz`
  and `cohort.vcf.gz` still overflow). The per-sample filtered VCFs keep the
  overflow until a re-call: the report bounds them on read, any other direct
  reader still sees it.
- **Report run time.** `own_span` + `_cohort_site_pass` take the dataset_paper
  family pages from 1 h 10 to 2 h 15 (2670 samples, 1206 families). Lead:
  resolve the site status of each own call once per sample (not per family),
  or join own calls and sites with a contig-sorted sweep rather than a
  per-call bisect.

### 1decies. Published evidence stores destroyed by publishDir (2026-09-29)

When dataset_paper was relaunched (filter -> report, large-DEL fix), the 2670
SVENGERS_FILTER tasks finished in 10 min with exit 0 — and **2071 of the 2670
published `samples/<S>/<S>.evidence.zarr` became unreadable**
(`OSError 40: Too many levels of symbolic links`).

Mechanism: `SVENGERS_FILTER` receives the store as input and **re-emits it as
output** (to chain it to genotype). Nextflow stages the input as a symbolic
link to the published path, then `publishDir` republishes *all* outputs —
including the store — to that same published path. Result: the published
path becomes a link to the work copy, which is itself a link to the published
path. A cycle. A directory cannot be hardlinked, so `mode: 'link'` silently
falls back to a symbolic link and nothing warns.

No data lost: the real stores live in the `SVENGERS_CALL` work directories,
only the pointer had broken. Repair: repoint each broken link to the real
store (chosen by number of contigs, then mtime), 2071/2071 repaired, 0
unrepairable. Checks: all 2670 stores open, the `sample` attribute matches
the directory for all 2670 (no sample mix-up), and a `svengers filter` re-run
on three repaired samples reproduces the published VCF **bit for bit** — the
repointed store is indeed the one the filter had read.

Fixes: `publishDir ... saveAs: { fn -> fn.endsWith('.evidence.zarr') ? null : fn }`
in `modules/svengers_filter.nf`, and `tests/test_workflow_modules.py`, which
fails if a module republishes an input that it re-emits (checked by restoring
the old version: the test fails).

Residual risk, not addressed: the published stores are links into the
`workflows/nextflow/work` tree (2.6 TB as of 2026-09-29). A `nextflow clean`
would destroy them, and they would have to be re-called from the CRAMs
(~6 days for 2670). Two options the day it matters: publish the stores in
`copy` mode (disk cost), or not publish them at all and treat the work dir as
the canonical location.

### 1nonies. Large deletions fell below the SQ threshold (2026-09-28)

Alexandre, on FAM-01-001: IGV shows a clean terminal deletion
chr22:41.9-50.8 Mb in the proband, SVengers only reported pieces of it.

The own calls did cover the whole region — but with SQs of 3.73 and 4.31 on
the two largest segments (6.43 Mb and 229 kb), hence `LowSQ`, hence absent
from the merge. Measured on 25 PMS genomes: **50% of 1 to 10 Mb DELs have
SQ < 5**.

Cause, in two steps:
1. The model has barely learned at these sizes: 131 examples > 1 Mb out of
   132,555 labeled DELs (0.10%), and **zero** true DUP examples above
   100 kb. Worse, these 131 have degenerate features (dhffc = 0,
   fc_chrom = 0, CN not called): they are zones of total coverage loss, not
   heterozygous deletions.
2. When the HMM splits a large deletion into adjacent segments, **the flanks
   of each segment fall within the neighboring segment**, so `dhffc` ≈ 1:
   the model's most discriminating feature goes blind. In the real case,
   dhffc = 1.03 whereas `fc_chrom` = 0.528 clearly says "one copy".

Fix: a DEL of at least `large_del_rescue_size` whose `fc_chrom` confirms the
loss (< 0.75) and that is lightly masked keeps its PASS regardless of SQ,
with `INFO/SQ_OVERRIDE=LARGE_DEL_DEPTH` for traceability. Threshold set on the
labeled data: at 200 kb the rescue takes 5 calls, all true; at 100 kb it
takes 23, of which 10 false; at 1 Mb it misses the first 229 kb segment of
the PMS terminal deletions. On SMP-07, 2 rescues over the whole genome, and
the PASS extent goes from scattered fragments to 41,894,001-50,799,001.

To do: this is a safeguard, not a fundamental fix. We need (a) realistic
examples of large events in v6 training — the PMS terminal deletions provide
some, provided they are validated independently; (b) to fix the flank
contamination of `dhffc` (exclude from the flanks the bins covered by another
call of the same sample); (c) nothing is done for large DUPs, for lack of a
single true example in training.

### 2bis. More discriminating impact score (Alexandre's request, 2026-09-28) — v1 DONE (2026-09-30)

**v1 in place** (`_GeneIndex.impact`, `_gene_points`, `_GeneTrek`): gene
content points following the ACMG/ClinGen CNV grid (Riggs et al. 2020,
sections 1-3), separate for DEL / DUP, with GeneTrek standing in for the
ClinGen dosage map. "Haploinsufficient-type" gene = established NDD
(`hc_ndd_v7` or `sysndd_definitive`) with dominant inheritance (any of the
three `inheritance*` columns). DEL: complete 1.00 / exonic 0.90; established
NDD with unknown mode 0.45; recessive only 0.15; candidate NDD
(candidate_ndd, SFARI 1, autism, PanelApp ID/epilepsy/autism) 0.30, 0.45 if
constrained; constrained only (`is_shet_constrained`, s_het ≥ 0.1 or
LOEUF < 0.35) 0.15; section 3: ≥ 25 / ≥ 35 coding genes +0.45 / +0.90. DUP:
intragenic DUP of a dominant NDD gene 0.45 (2I, possible loss of function),
complete DUP 0.15 at most (triplosensitivity absent from GeneTrek); ≥ 35 /
≥ 50 genes +0.45 / +0.90. No coding gene −0.60 (1B). Label: points, level,
rule, driver gene with its evidence (established/candidate NDD, AD/AR,
s_het), other NDD genes, number of coding genes. Numeric sort in the family
pages (`data-impact`).

Measured on the 576,429 PASS rows of the PMS family pages: the old score had
5 values and **3,664 rows at the maximum level**; the new one has 10, and
**49 rows at +1.90 — exactly the terminal 22q13 deletions removing SHANK3
with ≥ 35 coding genes**. Then come +1.45 / +1.00 (shorter SHANK3
deletions), +1.05 (8q24.3 duplications of 75-103 genes from the unbalanced
translocations, FAM-09 and FAM-06), then 499 exonic hits of dominant NDD
genes (+0.90, many without an own call).

Remaining for a v2: triplosensitivity (ClinGen dosage map, to be added to the
file), 5'/3' of the MANE transcript for partial DELs (2C vs 2D-2E), and
sections 4-5 (de novo, segregation, cohort phenotype) — to be calibrated on
the de novo already validated by IGV.

**Initial wording of the request:**


The current score is coarse: 0 intergenic / 1 intronic / 2 exonic, +1 if
`candidate_ndd_v5` (7316 genes), +2 if `hc_ndd_v7` (2200) — i.e. 4 possible
values and 2 columns read out of the **71** of
`genetrek-data-2026-09-15.tsv`. At this granularity, half of the coding
genome ends up at the same level.

Columns available and directly usable (fill rate measured on the 45,045 genes
of the file):

- **constraint / dosage**: `loeuf` (17,961), `shet_mean` (16,486),
  `is_shet_constrained` (3946), `constraint` (4475), `moeuf`, `soeuf`.
  This is the right prior for a **DEL** (loss-of-function intolerance).
  Caution: LOEUF/s_het say nothing about **triplosensitivity** — a DUP must
  not be scored on the same scale (lead: ClinGen triplosensitivity, absent
  from the file, to be added).
- **phenotype panels**: `panelapp_eng_id_*` (1514), `panelapp_aus_autism_*`,
  `*_epilepsy_*`, `sysndd_definitive_*`, `spark_*`, `sfari_1_*` (244),
  `autism_combined_v7`, `acmg_sf_v3_3` (secondary findings to report).
  Allows prioritization **according to the cohort** instead of a single
  score.
- **mechanism**: `chromatin_remodeling`, `epifactors_v2_1`, `syngo_v1_3`,
  `trrust_v2_tf` / `rummagene_tf`, `chrom_epi_tf_combined`.
- **mode of inheritance**: `inheritance`, `inheritance_panelapp_*`,
  `inheritance_sysndd_*` — an AD gene hit by a heterozygous DEL does not
  weigh the same as an AR gene.
- **quantitative autism evidence**: `37365347_autism_or`, `_p_value`,
  `_number_of_carriers_among_autistic_individuals`.
- `mane_select` (19,314) to decide cleanly on exonic involvement, and
  `omim_id` for the clinical link.

Proposal: a continuous score, separate by type (DEL/DUP), combining
(a) the fraction of the MANE transcript hit, (b) the constraint of the most
constrained overlapped gene, (c) membership in the panel of the cohort's
phenotype, (d) the mode of inheritance — rather than an integer from 0 to 4.
Keep a readable label next to the number, as today. To be calibrated on the
de novo lists already validated.

### 2ter. Genotype-level verification — DONE, but not duphold-style (2026-09-28)

Tested before implementing, on the 77 PMS trios, with the apparent de novo
rate as a proxy for genotyping error. **Two duphold approaches were measured
and discarded**:

- **DHFFC (fold-change against the flanks)**, the core of duphold: it
  discriminates nothing here, it is even inverted for DELs (1.4% error when
  it is "consistent" against 0.4% otherwise). The reason: **98% of DUP
  carriers are at sites smaller than 3 kb**, where depth has no resolution
  (1 kb bins) — the fold-change there is noise.
- **Requiring clips on both sides** of the breakpoint: inverted too (DUP
  10.5% error with both sides against 5.8% with one only).

What works, measured:

| criterion (among carriers GQ>=20) | DEL kept / error | DUP kept / error |
|---|---|---|
| none | 100% / 0.8% | 100% / 1.9% |
| **CN consistent with type** | 88% / **0.3%** | 24% / **0.3%** |
| RD consistent | 92% / 0.6% | 26% / 0.6% |
| CN consistent or PE>=2 | 92% / 0.4% | 68% / 1.3% |

The **re-genotyped CN** beats RD and costs nothing (already loaded by the
report). The junction channel is confirmed to be poorly specific: "SR only"
is the worst category (small DELs 5.9% error against 0.8% with PE>=2; small
DUPs 6.9% against 3.0%).

Done: **support** column ("strong" = a carrier's CN confirms the type) on the
events pages, with a selector. **Not** a default filter: on the 42 sites
reviewed in IGV, requiring it drops 2 of the 3 true ones (precision 23%
-> 20%). The two measurements contradict each other on this point — the trio
proxy says it is a good indicator of genotyping error, the IGV sample (3 true
ones, hence very fragile) says it discards true calls. To be settled when
there are more validations.

Still to do: make `_junction_support` specific (pair the clips at the
junction — the store does not keep the SA partner, so this requires adding it
to extraction), and a DHBFC (bins of comparable GC) if depth ever becomes the
limiting channel again on large events.

### 2quater. (archive) initial wording of the duphold request

What duphold (brentp) does: for each call, a depth fold-change against the
flanks (DHFFC), against the genome median (DHFC) and against bins of
comparable GC (DHBFC), plus junction support.

**We already do it, but in the wrong place.** `filter/features.py` computes
`dhffc` (exactly DHFFC), `fc_chrom` (~DHFC), `fc_residual`, `depth_cv`,
`frac_consistent`, and the model uses them — but **only on the per-sample
call**. Cohort re-genotyping, on the other hand, does no check of this kind:
this is precisely where the errors found on 2026-09-28 occur (96% of DUP
carriers judged false in IGV have RD = 1.00, strictly normal depth).

To do:
- move the check into `cohort/genotyper.py`: reject a carrier genotype whose
  depth contradicts the type when the size makes depth informative
  (measured: apparent de novo 10.1% for DELs with inconsistent depth against
  1.1% for consistent ones, 5.3% vs 2.1% for DUPs);
- add the **DHBFC** (comparison to bins of comparable GC), the only one of
  duphold's three that we lack — useful where GC correction leaves a residual
  bias;
- expose these values per genotype in the cohort VCF (FORMAT) and in the
  events pages, so that they are readable during IGV review;
- at the same time fix `_junction_support`, which counts clips near the
  breakpoint without checking that they support *this* junction: in the
  validation, all false carriers had junction support, so this channel
  discriminates nothing as it stands.

The GQ floor set on 2026-09-28 (`--min-carrier-gq`) treats the most visible
symptom; these three points treat the cause.
### 1octies. Cohort genotyper specificity: GQ floor (2026-09-28)

Large-scale measurement on PMS (77 trios, 32,386 sampled sites) of the
**apparent de novo rate** in the carrier child when both parents are
confident 0/0 — a direct proxy for the genotyping error rate:

| type | subgroup | carriers | apparent de novo |
|---|---|---|---|
| DEL | GQ < 20 | 20,832 | **16.7%** |
| DEL | GQ >= 20 | 251,294 | **0.8%** |
| DUP | GQ < 20 | 71,508 | **13.9%** |
| DUP | GQ >= 20 | 253,823 | **1.9%** |
| DEL | inconsistent depth | 27,069 | 10.1% |
| DUP | inconsistent depth | 254,989 | 5.3% |

GQ dominates every other criterion tested (depth consistency, size,
junction support — the latter discriminates nothing, all the false ones have
it, a sign that split-read counting is too permissive: it counts clips near
the breakpoint without checking that they support *this* junction).

Done: on the events pages, "carrier" = alternative genotype **and**
`GQ >= --min-carrier-gq` (default 20). Less certain genotypes remain
displayed in the per-member columns but no longer count as carrier status.

Check against the 42 sites reviewed in IGV (family FAM-04):

| criterion | true | false | precision |
|---|---|---|---|
| raw | 3/3 | 39/39 | 7% |
| GQ>=20 | 3/3 | 18/39 | 14% |
| backed by an own call | 3/3 | 18/39 | 14% |
| **GQ>=20 AND backed** | **3/3** | **10/39** | **23%** |

No true one lost, 29 of 39 false ones removed. Still to handle:
- the specificity of SR/PE counting in `_junction_support` (match the
  clips to the junction, not only to the position);
- the cohort statistics (NCARRIER/AF/CARRIER_FRAC) still count all
  alternative genotypes, with no GQ floor — to be decided, since touching
  them shifts AF, HWE and the derived tags.

### 1septies. Carrier status from re-genotyping alone is almost always false (2026-09-28)

**Blind** IGV validation by Alexandre on 60 events of family FAM-04 in PMS
(40 sites released by the directional criterion + 20 controls already PASS,
shuffled, FILTER hidden), 42 reviewed:

| group | n | true | false |
|---|---|---|---|
| released | 25 | 3 | 22 |
| PASS controls | 17 | 0 | 17 |

The controls at 0/17 say that the problem predates the directional change.
The cause lies elsewhere, and it is clear-cut — cross-referencing with the
family members' **own** calls:

| verdict | own call PASS | own call filtered | no own call |
|---|---|---|---|
| true (3) | 3 | 0 | 0 |
| false (39) | 6 | 12 | 21 |

A cohort site is a group of breakpoints; re-genotyping then assigns a
genotype to **all** samples, whether or not they produced a call. On family
FAM-04 in PMS: **80.4% of the 19,416 carried PASS events have no own call
behind them** (11.8% have a PASS call, 7.8% a filtered call). The review
scores 0/33 for the unbacked ones versus 3/9 for those backed by a PASS
call.

This explains the gap with the HG002 benchmark (precision 0.95 at SQ>=5):
the latter measures **per-sample calls**, not the re-genotyped carrier
status displayed in the events pages. And it explains why the same family
looked good in its standalone run (9 samples: few sites from others to
inherit).

Done: "own call" column (PASS / filtered / none) and selector on the events
pages, set by default to "with own call" — 19,416 -> 3,806 rows on family
FAM-04. The default keeps filtered calls, otherwise the rescued sites (the
familial SHANK3 DEL, LowSQ in every carrier) would be hidden.

To do: look at the specificity of the cohort genotyper itself (the
threshold for calling a carrier at re-genotyping) — that is the source of
the problem, the column only makes it visible.

### 1sexies. Genotyping discordance made directional (2026-09-28)

On the PMS cohort (Phelan-McDermid, 319 genomes), the causal SHANK3
deletion — chr22:50,654,001, 33 kb, the one in family FAM-04 — came out
`CohortArtifact`. Neither frequency nor HWE (p=1.0): it is `GTR_DISC`. Only
6 samples had nominated it at calling (low SQ, rescue); re-genotyping finds
64 carriers out of 319, relative gap 9.7 > threshold 0.5.

Yet this is the intended behavior of the pipeline: re-genotyping from the
stores is *meant* to be more sensitive than per-sample calling. The
criterion penalized this gain. Asymmetry measured on PMS: 36,328 sites
filtered in the "re-genotyping finds more" direction versus 425 in the
other.

The criterion is now **directional**: we only flag
`(nominations - porteurs) / nominations > seuil`, i.e. calls that
re-genotyping does not reproduce. The opposite artifact case ("the
genotyper calls everyone") has been covered since 2026-09-24 by
`HighCarrierFrac`, which did not exist when `GTR_DISC` was written.

Effect on PMS: 7,486 sites become PASS (+16%). Mendel comparison
(conditional on the child being a carrier, see below): MVR 22.9% for the
released sites versus 22.9% for the current PASS and 3.0% for the other
artifacts — they behave like the accepted set, not like the artifact class
(92% of which are carried by everyone). Caveat: their apparent de novo rate
is higher (14.5% versus 9.7%). Blind IGV validation in progress to decide.

**Reading trap in `mendel_report.tsv`**: the file only contains the
"informative" rows (`mendel.py`, "only report informative rows"), so
non-carrier children are absent from it. Any transmission rate computed
directly from this TSV is ~1.0 by construction — use the columns of the
summary or `bench trios`, which count over all genotyped sites.

### 1quinquies. `lcr_ovl` had been dead from the start (2026-09-25)

Alexandre: "events in an LCR are sometimes more often false, especially
DUPs". Check: `lcr_ovl > 0` for **0 calls out of 38,094** (whole genome)
and **0 out of 29,056** (HG002, training data). Cause: the
`selfchain-LCR.hg38.bed` track uses Ensembl naming (`1`, `2`…) whereas the
calls are in UCSC naming (`chr1`) — no overlap could ever match.
`IntervalSet` now resolves both conventions (test in test_reference.py).

Consequences:
- **No impact on v5**: since the feature was constant at 0 during training,
  no tree splits on it — verified, SQ strictly identical on 5000 calls with
  a random `lcr_ovl` (max difference 0.0).
- **For v6**: the feature becomes genuinely informative; retrain for it to
  be used. It is also a lesson about contracts: a constant feature should
  fail loudly — **safeguard added** on 2026-09-25:
  `train.check_feature_variance` refuses to train on a zero-variance (or
  entirely NaN) feature, with `--allow-constant-features` to force it.
  Per-type exemptions: `is_dup` (one model per type) and, on DUP, the
  DEL-specific channels (`cig`, `ins_z`, `ins_pairs`). Verified on the v5
  data: the safeguard would have refused training, naming `lcr_ovl`. The
  constant features detected are also recorded in the model manifest.

The LCR × DUP hypothesis itself: **not confirmed** on the 45 IGV
validations available once the feature was repaired — false DUPs with
`lcr>0` 3/8 (median 0.000) versus true DUPs 5/10 (median 0.008). What
really separates them, on this small sample: the type (DEL 26 true/1 false
= 96% vs DUP 10/8 = 56%) and `segdup_ovl` (3/8 of the false DUPs at 1.00,
0/10 of the true ones).

### 1quater. Per-sample genotyping safeguard (2026-09-25)

Follow-up to 1ter: recomputing the de novo safeguard while excluding
`HighCarrierFrac`/`CNDiscordant` changed nothing (already mutually
exclusive with PASS — 0 sites combine the two). The real cause of the
largest outliers (SMP-05/U/X, 12k-14k candidates each): three **unrelated**
trios with globally unreliable genotyping (60-69% of PASS carriers at
GQ&lt;20, median 10-13) — a sample problem, not a site problem.

Two safeguards added:
- `qc.WAVE_WARN` lowered from 0.145 to 0.128 (Youden) — captures SMP-06/X
  (0.139/0.140), below the old threshold.
- New (corrected twice on 2026-09-25): **carrier load**
  `report.carrier_excess` = number of **PASS** sites carried / cohort
  median, threshold **5x**, applied only from 20 samples upward (relative
  criterion). dataset_paper distribution: median 1.00, p90 4.4x, p95 6.4x,
  p99 9.6x, max 23.5x. The 5x threshold flags 212/2670 (7.9%) and captures
  18 of the 20 worst de novo outliers (median load of the outliers 1.41x
  versus 0.91x for the others).

  **Second calibration error** (same class as the first): the 1.3x
  threshold had been calibrated on the number of carriers over **all
  sites** (median 308,177) whereas the safeguard measures **PASS** carriers
  (median 3,284) — two unrelated distributions, hence 1037/2670 flagged on
  the first attempt. Rule to keep: calibrate on *exactly* the quantity the
  code computes, by reading the column produced by the code, never a
  neighboring quantity.

  **First, erroneous version**: the fraction of PASS carriers at GQ&lt;20,
  threshold 0.40 — calibration at real scale: 2670/2670 flagged. The cohort
  median of this fraction is 0.75 (p5 0.68, max 0.89) and the three real
  culprits (SMP-05/U/X) are at 0.61-0.69, i.e. percentiles 2-6: *among the
  cleanest*. What sets them apart is the load (487-497k carriers versus
  308k at the median, percentile 100). Lesson: do not calibrate a threshold
  on 3 cases without the background distribution. The GQ fraction is still
  displayed (it describes the callset: 75% of PASS carrier genotypes are
  below GQ 20) but is no longer used as a criterion.

Note: the hypothesis "it was `CommonRef` that inflated the metric on family
FAM-04" was wrong — the fraction is high everywhere, including on the
reference cohort without `--known-sites`. The current relative criterion,
on the other hand, remains inapplicable to small family runs (everyone
carries the same sites), hence the minimum of 20 samples.

### 1ter. v6 — clean the labels by transmission (finding of 2026-09-24)

IGV validation by Alexandre on family FAM-04 (v5, PASS sites): 10 false /
3 true among the first 13. The 10 false ones are carried by 8-9 members out
of 9 and are all `CohortArtifact` in dataset_paper (46-100% carriers out of
2670); the false DUPs have CN=1 in 98% of carriers.

Consequence for training: an artifact shared parent→child is labeled
**TP by transmission**. Test on the 16 wavy-genome trios: the rule "DUP
with dhffc<0.7" removes 22.6% of the "TP" and 26% of the FP —
indistinguishable, whereas on GIAB it touches almost nothing. The v5 labels
from wavy genomes are therefore contaminated by common familial artifacts,
and v5 learned to accept them.

To do for v6:
- `svengers label-trio --known-sites <cohort.vcf.gz>`: exclude from the
  labels (or turn into FP) the calls whose reference site has
  REF_CARRIER_FRAC ≥ 0.5, and those flagged HighCarrierFrac / CNDiscordant
  in the run;
- retrain, compare on HG002 (non-regression) and on the IGV validations
  (`validation` TSVs exported from the events pages — cross-reference them
  by site_id or by locus);
- on the post-processing side (done, 2026-09-24): `cohort-filter` sets
  HighCarrierFrac / CNDiscordant / CommonRef (`--known-sites`) and the merge
  propagates SQMAX/SQMEAN; the family run must point `--known_sites` at the
  latest large cohort.
- naming: `CohortArtifact` (HWE) also catches real common CNVs (site_199
  DEL 10.5 kb at 32% carriers, validated true) — consider a separate
  `CommonCNV` for sites in HWE disequilibrium but CN-consistent.

### 1bis. Detail of the initial training plan (reference)
- Download (proxy configured) into `resources/giab/`:
  - HG002/HG003/HG004 CRAMs ~35x GRCh38 (GIAB FTP, PCR-free ILMN)
  - HG002-Q100 v1.1 truth (`.../T2T-Q100_SVs_v1.1/GRCh38/`) + CMRG v1.00 SV
  - GIAB genome-stratifications (segdup/LCR/TR) + UCSC simpleRepeat + Umap
    k100
- Permissive `svengers call` on HG002/3/4 → labels by `truvari bench`
  (`--pctseq 0 --pctsize 0.7 --refdist 1000`, restricted to the confidence
  BED; outside the BED = excluded from training).
- 2 separate HistGradientBoosting models DEL / DUP, monotonic constraints
  (P(true DEL) decreasing in DHFFC, increasing in SR), CV by chromosome
  groups, out-of-fold isotonic calibration, SHAP → feature pruning.
- Delivery: `models/v1/` (joblib ×2 + manifest.json {version, schema_hash,
  data URLs, threshold→FDR table}). Target: FDR ≤ 2.8% DEL / 7.3% DUP
  (= gnomAD-SV v4).
- GHFC trio augmentation: NOT as v1 labels (circularity); MVR/transmission
  rate = evaluation metrics only.

### 2. Benchmark harness (`bench/`) — produces the paper's figures
- `bench truth`: Truvari bench+refine vs Q100/CMRG, size strata
  {50-300, 300-1k, 1-10k, 10-100k, 100k-1M, >1M} × context {clean, segdup,
  LCR, TR}; TSV+JSON output.
- `bench compare`: Manta/DELLY/dysgu (installed in tools/) + existing
  WisecondorX/ClinSV, same protocol; SVengers before/after filter.
- `bench trios`: MVR and transmission as a function of the SQ threshold on
  the GHFC cohorts (specificity figure independent of GIAB).
- `bench cohort` (2026-11-09): deployment gate for a model. Full cohort
  re-filtered and re-genotyped with the candidate, then PASS de novo
  GQ >= 20 per child (type × size) and recurrent de novo sites against
  production; alarm if the median rises by 25% and 2 calls/child.
  Rationale: v12a (v5r2 + 129 IGV verdicts ×20) had a per-genome panel at
  DUP -9.7% and clean "parental call" trios, but in the cohort DUP
  ≥ 1 kb de novo per child went from a median of 1 → 12 (GC coverage waves
  shared by the family). Mandatory criterion before any model change.
- **Platinum Pedigree (CEPH-1463), to do**: truth independent of GIAB for
  the paper. Reference: Genome Biology 2026,
  doi:10.1186/s13059-026-04219-3 (14 short-read callers on NA12878 and
  three of her children; DRAGEN F1 0.51, Manta 0.41, Dysgu 0.37; Truvari
  4.2.2 `--typeignore --passonly --pctseq 0`). HiFi + ONT truth for the
  family (`merged_hg38.svs.sort.oa.vcf.gz`, Platinum-Pedigree-Datasets
  repository): 13,460 DEL, 22,146 INS, 56 INV ≥ 50 bp, high-confidence BED
  for NA12878 only. Value: our models are trained on HG002/3/4, so an F1
  outside GIAB is more defensible; a pedigree with long-read truth to
  calibrate transmission and de novo. Plan:
  1. download the Illumina BAMs and the truth (Alexandre's approval
     required, sizes announced beforehand);
  2. produce the tracks in GHFC format (mosdepth 1 kb, DeepVariant VAF, see
     the track recipe) then `svengers call` without PoN, and with a PoN
     built on the family members (the GHFC PoN does not hold for another
     centre);
  3. compare DEL alone by size + tandem DUP extracted from the INS
     (`bench dup-truth`), not the all-type F1 (62% INS, which SVengers does
     not call); Manta and Dysgu re-run identically (tools/);
  4. parent → child transmission of SVengers calls against the truth;
  5. limitation to write up: the truth is rich in 100-500 bp, poor in
     > 10 kb, segdups and multi-allelic sites — it complements HG002, it
     does not replace it.

### 3. Hidden split reads — DONE (2026-09-08)
- `breakpoints/resolve.py`: fitting realignment (numba, no dependency) of
  the stored clips (store schema 0.2), clip→junction geometry, HOMLEN
  micro-homology, sub-bin dilution in the genotyper. Validated on chr20:
  exact ClinSV breakpoints + previously unseen small DUPs of 50-105 bp.
- Remaining for v1.1: recovery of mismatch-rich ends (NM/MD) without a
  clip, and adding HOMLEN as a feature of the ML filter (schema v2).

### 4. Panel of normals — DONE (2026-09-08)
- `cohort/pon.py`: streaming per-bin median+MAD (flat memory at 10k
  samples), sex-stratified, `svengers call --pon`, artifact catalogue
  (`pon artifacts` + `filter --artifact-bed`).
- Semantics validated on GHFC chr20 (panel of 20 EAGER founders): depth
  segments 264 → 14, 12/14 survivors rare. COMMON CNVs drop out of the
  depth channel (median/MAD) and remain detectable only through junctions
  → --pon reserved for rare-variant analyses, not for AF catalogues.
- Remaining for v1.1: optional SVD (k≈20) à la SavvyCNV; whole-genome PoN
  from the ~2000 GHFC founders (SLURM array over the stores).

### 5. Industrialization
- Nextflow wrapper modelled on tools/ghfc-ngs (DSL2, slurm+apptainer
  profiles).
- ~~INV then BND/translocations~~: `svengers bnd` (2026-11-06, see above).
- Packaging: PyPI (`svengers` — reservation to do, ALEXANDRE'S APPROVAL
  required before any publication), bioconda, container, mkdocs,
  CITATION.cff, Zenodo.
- Perf: profile `extract/scan.py`; if >40% of the runtime, a targeted Rust
  port (PyO3/rust-htslib) of this module alone.

## Locked design decisions (do not re-open without a reason)
- Cohort merge = bounded-centroid clustering, NO transitive chaining
  (SURVIVOR lesson: AF inflation 1.7-3.6×).
- Cohort re-genotyping from the stores ONLY (UKB svimmer+GraphTyper2
  pattern); never full joint calling.
- Filter = calibrated, interpretable gradient boosting, NO deep learning
  (bioRxiv 2026 Sedlazeck: simple RF = ResNet/Evo2 at F1 95.7%).
- Unified RD+PE/SR+BAF genotype likelihood = the paper's no. 1
  methodological claim; BAF also in the HMM emissions.
