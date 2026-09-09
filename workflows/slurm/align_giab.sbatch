#!/bin/bash
#SBATCH --job-name=giab-align
#SBATCH --partition=ghfc
#SBATCH --qos=ghfc
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/giab-align.%j.log
#
# GHFC-identical alignment of a GIAB sample's 40x NovaSeq PCR-free fastqs.
#   sbatch workflows/slurm/align_giab.sbatch HG003

set -euo pipefail
SAMPLE=${1:?usage: align_giab.sbatch SAMPLE}

ROOT=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
REF=/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta
FQ=$ROOT/resources/giab/fastq
OUT=$ROOT/resources/giab/${SAMPLE,,}
BWA=/pasteur/helix/projects/ghfc_wgs/tools/bwa-mem2-2.2.1/bwa-mem2
SAMBLASTER=/pasteur/helix/projects/ghfc_wgs/tools/bin/samblaster
module load samtools/1.22.1 2>/dev/null || export PATH=/opt/gensoft/exe/samtools/1.22.1/bin:$PATH

mkdir -p "$OUT"
TMP=${TMPDIR:-/local/scratch}/giab-align.$SLURM_JOB_ID
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS=$SLURM_CPUS_PER_TASK

$BWA mem -Y -t $((THREADS - 8)) \
    -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA\tLB:novaseq-pcr-free-40x" \
    "$REF" \
    "$FQ/$SAMPLE.novaseq.pcr-free.40x.R1.fastq.gz" \
    "$FQ/$SAMPLE.novaseq.pcr-free.40x.R2.fastq.gz" \
  | $SAMBLASTER --addMateTags --ignoreUnmated \
  | samtools sort -@ 8 -m 3G -T "$TMP/sort" -O bam -o "$TMP/$SAMPLE.sorted.bam"

samtools view -@ $THREADS -C -T "$REF" -o "$OUT/$SAMPLE.GRCh38_GIABv3.cram" "$TMP/$SAMPLE.sorted.bam"
samtools index -@ $THREADS "$OUT/$SAMPLE.GRCh38_GIABv3.cram"
samtools flagstat -@ $THREADS "$OUT/$SAMPLE.GRCh38_GIABv3.cram" > "$OUT/$SAMPLE.flagstat.txt"
echo "DONE $(date)"
