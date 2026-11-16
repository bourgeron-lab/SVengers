#!/bin/bash
#SBATCH --job-name=giab-call
#SBATCH --partition=ghfc
#SBATCH --qos=ghfc
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=logs/giab-call.%j.log
# svengers call + filter (features) for one GIAB sample.
#   sbatch workflows/slurm/call_giab.sbatch HG003
set -euo pipefail
SAMPLE=${1:?usage: call_giab.sbatch SAMPLE}
ROOT=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
V=$ROOT/.mamba/envs/vedette/bin/vedette
OUT=$ROOT/resources/giab/${SAMPLE,,}_run
mkdir -p "$OUT"
$V call --sample $SAMPLE --out-dir "$OUT" --config $ROOT/configs/giab.toml \
    --resources $ROOT/resources \
    --cram $ROOT/resources/giab/${SAMPLE,,}/$SAMPLE.GRCh38_GIABv3.cram
$V filter "$OUT/$SAMPLE.svengers.vcf.gz" \
    --store "$OUT/$SAMPLE.evidence.zarr" \
    --out "$OUT/$SAMPLE.tier1.vcf.gz" \
    --config $ROOT/configs/giab.toml \
    --features-out "$OUT/$SAMPLE.features.tsv"
echo "DONE $(date)"
