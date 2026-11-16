# Running SVengers on the cohort via SLURM

The cohort pipeline is sample-parallel by design. Typical GHFC run:

```bash
SVENGERS_DIR=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
COHORT=EAGER   # any cohort under WGS/GHFC-GRCh38/cohorts/

# 1. manifest from the pedigree
$SVENGERS_DIR/.mamba/envs/vedette/bin/vedette manifest discover \
    /pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples \
    --pedigree /pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts/$COHORT/$COHORT.pedigree.tsv \
    --out $COHORT.manifest.tsv

# 2. per-sample call+filter as a SLURM array (one task per sample)
sbatch --array=1-$(( $(wc -l < $COHORT.manifest.tsv) - 1 ))%200 \
    $SVENGERS_DIR/workflows/slurm/call_filter_array.sbatch $COHORT.manifest.tsv out/$COHORT

# 3. merge, per-sample re-genotyping (second array), assembly + pedigree
#    (see call_filter_array.sbatch header for the follow-up commands)
```

Notes
- Steps are idempotent: re-submitting skips nothing but overwrites atomically;
  restrict with `--array` ranges to retry failures.
- `svengers genotype` reads only the evidence stores (~100 MB/sample), so the
  re-genotyping array is I/O-light and can run at high concurrency.
