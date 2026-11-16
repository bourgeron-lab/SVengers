#!/bin/bash
#SBATCH --job-name=svengers-call
#SBATCH --partition=common
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --output=logs/svengers-call.%A_%a.log
#
# Per-sample SVengers call + filter, one SLURM array task per manifest row.
#
#   sbatch --array=1-<N>%200 call_filter_array.sbatch <manifest.tsv> <out_dir>
#
# Follow-up (after the array finishes):
#   svengers merge out/*/**.filtered.vcf.gz --out sites.vcf.gz --fasta $FASTA
#   sbatch --array=1-<N> genotype array (svengers genotype sites.vcf.gz --store ...)
#   svengers cohort-filter sites.vcf.gz --genotypes @gt_list.txt --out cohort.vcf.gz --pedigree $PED
#   svengers ped-annotate cohort.vcf.gz --pedigree $PED --out-report mendel.tsv

set -euo pipefail

MANIFEST=$1
OUT_ROOT=$2
SVENGERS_DIR=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
SVENGERS=$SVENGERS_DIR/.mamba/envs/vedette/bin/vedette
CONFIG=$SVENGERS_DIR/configs/ghfc.toml
RESOURCES=$SVENGERS_DIR/resources   # shared cached reference tracks

mkdir -p logs "$OUT_ROOT"

# manifest row for this array task (line 1 is the header)
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
SAMPLE=$(cut -f1 <<<"$ROW")
CRAM=$(cut -f2 <<<"$ROW")
DEPTH_BG=$(cut -f5 <<<"$ROW")
VAF_BG=$(cut -f6 <<<"$ROW")

OUT="$OUT_ROOT/$SAMPLE"
if [[ -s "$OUT/$SAMPLE.filtered.vcf.gz" ]]; then
    echo "$SAMPLE already done"; exit 0
fi

ARGS=(--sample "$SAMPLE" --out-dir "$OUT" --resources "$RESOURCES" --config "$CONFIG" --cram "$CRAM")
[[ -n "$DEPTH_BG" ]] && ARGS+=(--depth-bedgraph "$DEPTH_BG")
[[ -n "$VAF_BG" ]] && ARGS+=(--vaf-bedgraph "$VAF_BG")

$SVENGERS call "${ARGS[@]}"

$SVENGERS filter "$OUT/$SAMPLE.svengers.vcf.gz" \
    --store "$OUT/$SAMPLE.evidence.zarr" \
    --out "$OUT/$SAMPLE.filtered.vcf.gz" \
    --config "$CONFIG" \
    --features-out "$OUT/$SAMPLE.features.tsv"
