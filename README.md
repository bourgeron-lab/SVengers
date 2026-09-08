# Vedette

**Vedette** is a multi-evidence germline CNV caller (DEL/DUP, later INV/TRA) for short-read
WGS, designed around one priority: **specificity**. It jointly models read depth, discordant
pairs, split reads and B-allele frequency in a single genotype likelihood, and ships a
tiered, interpretable post-calling filter (hard filters + calibrated gradient-boosting score
trained on GIAB HG002-Q100).

Scales from 1 to 10,000 samples: per-sample calling writes a compact *evidence store*
(~100 MB/sample, Zarr) from which cohort-wide re-genotyping runs without ever re-reading a
CRAM. Native CRAM input, pedigree-aware (Mendelian consistency, de novo scoring).

## Subcommands

| command | role |
|---|---|
| `vedette call` | single-sample calling: CRAM → evidence store + VCF |
| `vedette filter` | tiered FP filtering: annotations + calibrated SQ score |
| `vedette merge` | cross-sample clustering of filtered calls → cohort sites |
| `vedette genotype` | force-genotype cohort sites in every sample from evidence stores |
| `vedette cohort-filter` | AF, HWE, artifact and batch metrics on the cohort VCF |
| `vedette ped` | Mendelian consistency, inheritance, de novo quality |
| `vedette pon` | panel-of-normals depth profile and artifact catalogue |
| `vedette manifest` | build/validate the cohort sample manifest |

Status: v0.1 in active development (Institut Pasteur, GHFC lab).

License: MIT.
