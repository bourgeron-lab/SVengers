#!/bin/bash
#SBATCH --job-name=hg002-align
#SBATCH --partition=common
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/hg002-align.%j.log
#
# Align HG002 NovaSeq PCR-free fastqs EXACTLY like the GHFC cohort
# (bwa-mem2 -Y + samblaster --addMateTags --ignoreUnmated, CRAM vs the GHFC
# reference), so the trained filter sees cohort-identical data.

set -euo pipefail

ROOT=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
REF=/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta
FQ=$ROOT/resources/giab/fastq
OUT=$ROOT/resources/giab/hg002
BWA=/pasteur/helix/projects/ghfc_wgs/tools/bwa-mem2-2.2.1/bwa-mem2
SAMBLASTER=/pasteur/helix/projects/ghfc_wgs/tools/bin/samblaster
module load samtools/1.22.1 2>/dev/null || export PATH=/opt/gensoft/exe/samtools/1.22.1/bin:$PATH

mkdir -p "$OUT"
TMP=${TMPDIR:-/local/scratch}/hg002-align.$SLURM_JOB_ID
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS=$SLURM_CPUS_PER_TASK

$BWA mem -Y -t $((THREADS - 8)) \
    -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tLB:novaseq-pcr-free-40x' \
    "$REF" \
    "$FQ/HG002.novaseq.pcr-free.40x.R1.fastq.gz" \
    "$FQ/HG002.novaseq.pcr-free.40x.R2.fastq.gz" \
  | $SAMBLASTER --addMateTags --ignoreUnmated \
  | samtools sort -@ 8 -m 3G -T "$TMP/sort" -O bam -o "$TMP/HG002.sorted.bam"

samtools view -@ $THREADS -C -T "$REF" -o "$OUT/HG002.GRCh38_GIABv3.cram" "$TMP/HG002.sorted.bam"
samtools index -@ $THREADS "$OUT/HG002.GRCh38_GIABv3.cram"
samtools flagstat -@ $THREADS "$OUT/HG002.GRCh38_GIABv3.cram" > "$OUT/HG002.flagstat.txt"
echo "DONE $(date)"
