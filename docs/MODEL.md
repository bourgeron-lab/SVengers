# Filter models v15 and v12b — model card

`models/v15` is the default filter model of SVengers since v2.3.0;
`models/v12b`, its predecessor, is the one frozen with v2.1.0 for the paper.
Both score each raw DEL/DUP call (SQ, a calibrated Phred-scaled probability
that the call is real) from the same evidence-store features
(`manifest.json`, feature-schema hash `f89b9e8890da996f`) with the same
recipe; v15 adds training rows from the cohort fast path (below). Two gradient-boosted classifiers (scikit-learn 1.7
`HistGradientBoostingClassifier`, monotone constraints where the physics is
unambiguous), one for DEL and one for DUP, each calibrated by isotonic
regression on out-of-fold probabilities.

## Training data

| Source | Access | DEL TP / FP | DUP TP / FP | Labels | Model |
|---|---|---|---|---|---|
| GIAB HG002, HG003, HG004 (NovaSeq PCR-free, about 40x), depth from the CRAM | public | 11,047 / 14,104 | 3,862 / 18,017 | HG002: Truvari against GIAB v5.0q; HG003/HG004 (its parents): by transmission against the HG002 truth set and HG002 false-positive loci (`svengers label --parent-mode`) | v12b, v15 |
| 16 children of GHFC family trios, chosen for wavy depth | controlled | 68,024 / 2,832 | 39,628 / 3,145 | transmission from both parents (`svengers label-trio`) | v12b, v15 |
| GIAB HG002, HG003, HG004 re-called on the cohort fast path (1 kb bedgraph depth, PoN variance mode, v2.2) | public | 10,776 / 7,921 | 3,609 / 16,874 | as the first row | v15 |

- The weights are public (`models/`); the 16 GHFC genomes are not. They
  come from research cohorts; access follows the cohorts' data-access
  procedure.
- The GHFC children supply 71 % of the labels and every example of a wavy
  genome (depth oscillating with GC content and replication timing). GIAB has
  no wavy genome.
- Transmission labels: a child's call is TP when a parent confidently
  carries the same event (genotype GQ >= 20), FP when both parents are
  confidently reference there.
- `wave_score` is pinned to -1 (neutral) at scoring: the model does not see
  the genome's waviness, it learns what a wave artefact looks like.

Out-of-fold metrics of v12b (all labels, 5 folds of whole chromosomes): DEL
PR-AUC 0.995 (n = 96,007), DUP 0.987 (n = 64,652); v15: DEL 0.993
(n = 114,704; at SQ >= 7, recall 0.914 and FDR 0.023), DUP 0.980
(n = 85,135; recall 0.824, FDR 0.035). PASS is SQ >= 5 (`[filter]
low_sq_threshold`); at SQ >= 7, DEL recall 0.936 and FDR 0.019, DUP recall
0.863 and FDR 0.029 (the manifest's threshold table). These are on the
training label sets and are not comparable across models trained on
different labels.

## Out-of-sample HG002

Retrained with the same recipe without HG002 (`v12b_noHG002`), on HG002:
DEL F1 0.515 against GIAB v5.0q (0.836 for 1-10 kb, precision 0.872) and
0.582 against CMRG.

Caveat: `v12b_noHG002` keeps the HG003/HG004 labels, which were defined from
the HG002 truth set and HG002 false-positive loci, so it is out of sample for
the genome, not fully for the loci. Locus-disjoint check (2026-11-17,
`runs/_chromsplit/`): two models trained without HG002, each on half of the
chromosomes (odd autosomes + chrX / even autosomes + chrY), HG002 scored on
the chromosomes its model never saw. DEL F1 against GIAB v5.0q: 0.515 overall
(precision 0.817, recall 0.376), 0.843 at 1-10 kb, 0.717 at 10-100 kb, 0.400
at 50-300 bp; CMRG 0.552; tandem duplications 0.130. The F1 is unchanged
overall, precision is 5 points lower; each half-model also had half the
labels, so this is a conservative estimate.

## HG002 in the cohort fast path

The numbers above use depth counted from the CRAM, with no panel of normals
(PoN). A cohort is usually run on the fast path: depth from the 1 kb mosdepth
bedgraph (`by1000`) normalized by a PoN built from cohort founders. Up to
v2.1, the PoN worked in **mask mode**: each bin divided by the founders' median
and masked where their depth varies across the population. Those bins are where
polymorphic CNVs lie: on HG002 the mask covers 8.2 % of the autosomal bins
(4.6 % without a PoN), but more than 70 % of the bins of 44 % of the true
deletions of 1-10 kb and of 65 % of those of 10-100 kb (none without it;
`runs/_cmp/hg002_pon_mask.py`).

Since v2.2 the GHFC configuration reads the PoN in **variance mode** (PoN
schema 0.2, `[segment] pon_mode = "variance"`). Per bin and sex stratum the
panel stores the centre and spread of the founders at the expected copy
number (the reference group), not their median: at a common deletion whose
non-carriers are a minority the median sits among the carriers, the reference
group does not. A bin is masked only when the founders have no reference group
or do not sit at integer copy numbers (a spread-out depth, not a CNV);
elsewhere the sample's ratio is divided by the group's centre and the group's
spread, in excess of the contig's typical spread (the founders' own sampling
noise, already in the model), enters the HMM and genotyper depth likelihoods
as extra variance. Bins where GRCh38 carries the minor allele (most founders
deleted, at integer states of the reference scale) are read against the
reference (`INFO/REFMINOR`).

HG002 re-called on the fast path, scored by `v12b_noHG002`, DEL F1 against
GIAB v5.0q:

| Path | 1-10 kb | 10-100 kb | CMRG |
|---|---|---|---|
| CRAM depth, no PoN (above) | 0.836 | 0.743 | 0.582 |
| bedgraph, no PoN | 0.839 | 0.718 | 0.567 |
| bedgraph + PoN mask mode, v2.1.0 | 0.465 | 0.370 | 0.455 |
| bedgraph + PoN mask mode, v2.1.4 | 0.565 | 0.370 | 0.480 |
| bedgraph + PoN variance mode, v2.2.0 | 0.705 | 0.706 | 0.496 |

The depth source does not matter; the PoN mode does. In mask mode, at
1-10 kb the PoN raises precision a little (0.934 without, 0.966-0.974 with)
and cuts recall from 0.76 to 0.31 (v2.1.0) or 0.40 (v2.1.4; v2.1.4 recovers
the deletions that the mask had hard-filtered as MaxN although split or
discordant reads support them). Most of that loss was common deletions: at
1-10 kb, recall is 0.86 without the PoN and 0.45 with it in mask mode for
deletions with gnomAD-SV allele frequency >= 5 % (n = 545), 0.47 and 0.14 for
deletions absent from gnomAD-SV (n = 116), 0.50 and 0.39 below 1 % (n = 46).
Variance mode brings 1-10 kb recall to 0.557 at precision 0.960, and
10-100 kb to 0.545 at precision 1.0 (CRAM path without a PoN: 0.739 and
0.591).

## Ablation: public data only (v14p)

Same recipe, features and labels, trained on GIAB HG002/3/4 alone
(25,151 DEL and 21,879 DUP labels).

- **HG002 out of sample** (trained without HG002): DEL F1 0.513 vs 0.515 for
  v12b (1-10 kb: 0.843 vs 0.836; precision 0.850 vs 0.872), CMRG 0.536 vs
  0.582. On a clean genome the public model is as good.
- **186 GHFC genomes**: PASS per genome DEL -2.8 %, DUP -25 % (mostly
  < 1 kb). In trios, the children lose about 520 DEL of 1-10 kb that v12b
  keeps, 97.9 % of them carried by a parent: real deletions.
- **2,670 GHFC genomes, wavy stratum** (62 genomes, wave score >= 0.20):
  deletions >= 10 kb just below PASS, the ones cohort recurrence can rescue,
  rise from 3.2 to 10.4 per genome. These are wave artefacts; rescued, they
  become false de novo calls. The pre-deployment gate fails.

Without wavy genomes in training, the model no longer rejects wave artefacts.
v14p is not shipped; it is the reason v12b keeps its controlled-access
training data.

## v15: training on the fast path too (v2.3.0)

The rows of the table above are v12b scoring fast-path calls it had never
seen: its training depth came from the CRAM (reads at MAPQ >= 1, every base),
the fast path reads the median depth of 1 kb bins at MAPQ >= 40. A deletion
of a few kb spans few informative bins there, its masked fraction and
depth-agreement features look like an artefact's, and v12b gave it a low SQ:
of the true deletions of 1-10 kb that the fast path called but filtered as
`LowSQ`, 149 / 149 reached SQ >= 5 once their features were replaced by the
same call's CRAM-path features (`runs/_loss_v22/counterfactual.py`).

v15 is v12b's training set plus GIAB HG002/3/4 re-called on the fast path
(third row of the table), same recipe. On HG002 out of sample (fast path,
v2.2 calls, GIAB v5.0q), against `v12b_noHG002`:

| DEL | v12b | v15 |
|---|---|---|
| all sizes, F1 (precision) | 0.461 (0.917) | 0.498 (0.900) |
| 50-300 bp, F1 (precision) | 0.346 (0.878) | 0.379 (0.842) |
| 300 bp-1 kb, F1 | 0.686 | 0.704 |
| 1-10 kb, F1 (recall, precision) | 0.705 (0.557, 0.960) | 0.820 (0.714, 0.964) |
| 10-100 kb, F1 | 0.706 | 0.706 |
| CMRG, F1 | 0.496 | 0.530 |

On the CRAM path the two are equal (DEL F1 0.513 vs 0.515, 1-10 kb 0.843 vs
0.836): v15 replaces v12b on both paths. A variant trained on the fast-path
GIAB rows without the CRAM-path ones lies on the same precision-recall curve
at a lower precision at SQ >= 5 and was not kept.

On the 2,989 GHFC genomes (v2.2.0 calls re-filtered), the cohort gate is
quiet in every wave stratum: median de novo per child 7 in both (6 in the 491
children of wave score < 0.128), recurrent de novo sites 495 -> 465,
confirmed de novo 182 -> 184 (dataset_paper) and 147 -> 147 (PMS). PASS per
genome (carrier genotypes, dataset_paper): DEL 1-10 kb 510 -> 532, DEL < 1 kb
1,875 -> 1,864, DUP < 1 kb 330 -> 302, larger events unchanged. Deletions and
duplications curated by eye in IGV: DEL 518 -> 519 / 562, DUP 306 / 410 and de
novo DEL 48 / 48 unchanged.

## Reproducing

The training recipe is `svengers train` on labeled feature TSVs
(`svengers label` for genomes with a truth set, `svengers label-trio` for
trio children), then `neutral_features = {"wave_score": -1.0}` in the
manifest. With other wavy trios (for example public high-coverage trios), the
same commands train an equivalent model; check it with `svengers bench
pregate` and `svengers bench cohort` before use.
