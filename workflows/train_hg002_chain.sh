#!/bin/bash
# Post-alignment training chain: HG002 CRAM -> calls -> features -> Truvari
# labels -> trained calibrated models -> filtered re-evaluation.
# Run on a compute node or locally (call step: ~15 min / 8 cores).

set -euo pipefail
ROOT=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
cd "$ROOT"
V=$ROOT/.mamba/envs/vedette/bin/vedette
GIAB=$ROOT/resources/giab
OUT=$GIAB/hg002_run
REF=/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta
TRUTH=$GIAB/HG002_GRCh38_v5.0q_stvar.vcf.gz
BED=$GIAB/HG002_GRCh38_v5.0q_stvar.benchmark.bed

mkdir -p "$OUT"

echo "== 1. call (self-contained CRAM scan, 1kb bins) =="
$V call --sample HG002 --out-dir "$OUT" --config configs/giab.toml \
    --resources "$ROOT/resources" \
    --cram "$GIAB/hg002/HG002.GRCh38_GIABv3.cram"

echo "== 2. filter (features only, no model) =="
$V filter "$OUT/HG002.vedette.vcf.gz" \
    --store "$OUT/HG002.evidence.zarr" \
    --out "$OUT/HG002.tier1.vcf.gz" \
    --config configs/giab.toml \
    --features-out "$OUT/HG002.features.tsv"

echo "== 3. label vs GIAB v5.0q stvar =="
$V label "$OUT/HG002.features.tsv" \
    --calls-vcf "$OUT/HG002.vedette.vcf.gz" \
    --truth-vcf "$TRUTH" --truth-bed "$BED" --fasta "$REF" \
    --truvari-dir "$OUT/truvari_label" \
    --out "$OUT/HG002.labeled.tsv"

echo "== 4. train models v1 =="
$V train "$OUT/HG002.labeled.tsv" --out models/v2 --model-version v2 \
    --training-data "HG002 40x NovaSeq PCR-free (brain-genomics-public), GHFC-style alignment, truth GIAB v5.0q stvar"

echo "== 5. re-filter with the trained model =="
$V filter "$OUT/HG002.vedette.vcf.gz" \
    --store "$OUT/HG002.evidence.zarr" \
    --out "$OUT/HG002.filtered_v2.vcf.gz" \
    --config configs/giab.toml \
    --model-dir models/v2

echo "== 6. stratified evaluation: tier1-PASS vs SQ-filtered =="
$V bench truth "$OUT/HG002.tier1.vcf.gz" \
    --truth "$TRUTH" --bed "$BED" --fasta "$REF" \
    --segdup-bed /pasteur/helix/projects/ghfc_wgs/references/GRCh38/clinsv/refdata-b38/tracks/GRCh38GenomicSuperDup.bed.gz \
    --lcr-bed /pasteur/helix/projects/ghfc_wgs/references/GRCh38/selfchain-LCR.hg38.bed \
    --out-dir "$OUT/eval_tier1" --plot || true

echo "DONE $(date)"
