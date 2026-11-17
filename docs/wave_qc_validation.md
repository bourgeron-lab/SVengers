# Wave QC — score validation (2026-09-14)

`svengers.qc.sample_wave_score` measures the fraction of low-frequency
spectral energy (periods >= 200 kb, 10 Mb windows, linear detrend) of the
**GC-corrected** coverage ratio — the signal that actually feeds the HMM.
External truth: classifications from the GHFC pipeline
`Wave_hg38_remapping_slurm_nextflow` (threshold = 95th calibration
percentile over 4017 raw bedgraphs).

## Protocol

592 genomes from the dataset_paper cohort with a published evidence store:
the 172 available genomes labelled WAVES + 420 NORMAL drawn at random
(seed 1). Score computed on `depth_corr / mean_depth` per contig (1 kb
bins, exclude/N masks applied), mean of the per-window medians.

## Results

| metric | value |
|---|---|
| AUC vs GHFC labels | **0.946** |
| Spearman vs GHFC composite score | 0.651 |
| median NORMAL / WAVES | 0.112 / 0.174 |
| p95 NORMAL | 0.145 |

Candidate thresholds (on this sample, enriched to 29% wavy genomes):

| threshold | sensitivity | specificity |
|---|---|---|
| 0.128 (Youden) | 0.93 | 0.87 |
| **0.145 (p95 of normals — retained, `qc.WAVE_WARN`)** | **0.80** | **0.95** |
| 0.160 | 0.67 | 0.97 |

The dynamic range is compressed compared with the scores of the GHFC
detector (computed on raw coverage): GC correction absorbs part of the
artefact. The operational link with runtime is documented in the same
analysis: WAVES samples have a median `svengers call` duration 5× longer
(8.5 h vs 1.7 h, odds ratio 8.8 for the slowest decile), hence the
`breakpoint.resolve_max_alignments` budget introduced at the same time as
this score.

Reproduce: `scratch validate_wave.py` (session), dataset_paper stores +
`final_results/merged_summary.csv` from the Wave pipeline.

## Update 2026-09-25: threshold lowered to 0.128 (Youden)

The de novo safeguard (report.py) let through two wavy-genome samples at
0.139 and 0.140 (SMP-06/X), below the p95 threshold of 0.145, with the same
profile as a confirmed case (SMP-05, 0.178) — globally unreliable
genotyping, thousands of false "de novo" calls in their child despite all
the site filters. The `qc.WAVE_WARN` threshold moves to the Youden point
already documented above (0.128; sensitivity 0.93, specificity 0.87), at
the cost of slightly more false positives on genuinely normal samples.
