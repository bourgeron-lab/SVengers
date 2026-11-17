#!/bin/bash
# v4 training chain: HG002/3/4 stores backfilled with wave_score (schema 0.4),
# re-extracted features (30-feature schema incl. wave_score), HG002 Truvari
# labels + HG003/4 transmission semi-supervision, HistGB v4, re-eval.
# LOGIN NODE ONLY: truvari SIGILLs on ghfc compute nodes.
set -euo pipefail
ROOT=/pasteur/helix/projects/ghfc_wgs/Divers/Alex/Software/CNV
cd "$ROOT"
V=$ROOT/.mamba/envs/vedette/bin/vedette
PY=$ROOT/.mamba/envs/vedette/bin/python
GIAB=$ROOT/resources/giab
REF=/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta
TRUTH=$GIAB/HG002_GRCh38_v5.0q_stvar.vcf.gz
BED=$GIAB/HG002_GRCh38_v5.0q_stvar.benchmark.bed
# raw calls of a run: <S>.svengers.vcf.gz, or <S>.vedette.vcf.gz before v2.0
raw_calls() { if [ -e "$1/$2.svengers.vcf.gz" ]; then echo "$1/$2.svengers.vcf.gz"; else echo "$1/$2.vedette.vcf.gz"; fi; }

echo "== 1. backfill wave_score dans les stores GIAB =="
$PY - <<'PYEOF'
import zarr
from svengers.extract.store import EvidenceStore
from svengers.qc import sample_wave_score

for run, s in (("hg002_run","HG002"), ("hg003_run","HG003"), ("hg004_run","HG004")):
    path = f"resources/giab/{run}/{s}.evidence.zarr"
    st = EvidenceStore(path)
    md = float(st.meta["mean_depth"])
    ratio, mask = {}, {}
    for c in st.contigs:
        _r, corr, m = st.depth(c)
        ratio[c] = corr / md
        mask[c] = m.astype(bool)
    score, by_contig = sample_wave_score(ratio, mask, st.bin_size)
    root = zarr.open_group(path, mode="r+")
    root.attrs["wave_score"] = score
    root.attrs["wave_by_contig"] = {k: round(v, 4) for k, v in by_contig.items()}
    root.attrs["schema_version"] = "0.4"
    print(f"{s}: wave_score={score:.4f}")
PYEOF

echo "== 2. features (schema 30 with wave_score) =="
for S in HG002 HG003 HG004; do
  run=$GIAB/$(echo $S | tr 'A-Z' 'a-z')_run
  $V filter "$(raw_calls "$run" "$S")" \
      --store "$run/$S.evidence.zarr" \
      --out "$run/$S.tier1_v4.vcf.gz" \
      --config configs/giab.toml \
      --features-out "$run/$S.features_v4.tsv" &
done
wait

echo "== 3. labels: HG002 truvari =="
$V label "$GIAB/hg002_run/HG002.features_v4.tsv" \
    --calls-vcf "$(raw_calls "$GIAB/hg002_run" HG002)" \
    --truth-vcf "$TRUTH" --truth-bed "$BED" --fasta "$REF" \
    --truvari-dir "$GIAB/hg002_run/truvari_label_v4" \
    --out "$GIAB/hg002_run/HG002.labeled_v4.tsv"

echo "== 4. labels: HG003/HG004 semi-supervision par transmission =="
for S in HG003 HG004; do
  run=$GIAB/$(echo $S | tr 'A-Z' 'a-z')_run
  $V label "$run/$S.features_v4.tsv" \
      --calls-vcf "$(raw_calls "$run" "$S")" \
      --truth-vcf "$TRUTH" --truth-bed "$BED" --fasta "$REF" \
      --truvari-dir "$run/truvari_label_v4" \
      --parent-mode \
      --recurrent-fp-vcf "$GIAB/hg002_run/truvari_label_v4/fp.vcf.gz" \
      --out "$run/$S.labeled_v4.tsv"
done

echo "== 5. train v4 =="
$V train \
    "$GIAB/hg002_run/HG002.labeled_v4.tsv" \
    "$GIAB/hg003_run/HG003.labeled_v4.tsv" \
    "$GIAB/hg004_run/HG004.labeled_v4.tsv" \
    --out models/v4 --model-version v4 \
    --training-data "HG002 (self-labeled) + HG003/HG004 (transmission semi-supervised), 40x NovaSeq GHFC-style, truth GIAB v5.0q stvar, 30-feature schema incl. wave_score"

echo "== 6. re-filter HG002 avec v4 + sous-ensembles PASS / SQ>=7 =="
$V filter "$(raw_calls "$GIAB/hg002_run" HG002)" \
    --store "$GIAB/hg002_run/HG002.evidence.zarr" \
    --out "$GIAB/hg002_run/HG002.filtered_v4.vcf.gz" \
    --config configs/giab.toml \
    --model-dir models/v4
$PY - <<'PYEOF'
import pysam
vin = pysam.VariantFile("resources/giab/hg002_run/HG002.filtered_v4.vcf.gz")
outs = {
    "pass": pysam.VariantFile("resources/giab/hg002_run/HG002.v4_pass.vcf.gz", "wz", header=vin.header),
    "sq7": pysam.VariantFile("resources/giab/hg002_run/HG002.v4_sq7.vcf.gz", "wz", header=vin.header),
}
for rec in vin:
    filt = set(rec.filter.keys()) - {"PASS"}
    sq = rec.samples[0].get("SQ")
    if not filt:
        outs["pass"].write(rec)
    if (not filt or filt == {"LowSQ"}) and sq is not None and sq >= 7:
        outs["sq7"].write(rec)
for o in outs.values():
    o.close()
for k in outs:
    pysam.tabix_index(f"resources/giab/hg002_run/HG002.v4_{k}.vcf.gz", preset="vcf", force=True)
PYEOF

echo "== 7. stratified eval v4 (SQ>=7 and PASS) =="
for SUB in sq7 pass; do
  $V bench truth "$GIAB/hg002_run/HG002.v4_$SUB.vcf.gz" \
      --truth "$TRUTH" --bed "$BED" --fasta "$REF" \
      --segdup-bed /pasteur/helix/projects/ghfc_wgs/references/GRCh38/clinsv/refdata-b38/tracks/GRCh38GenomicSuperDup.bed.gz \
      --lcr-bed /pasteur/helix/projects/ghfc_wgs/references/GRCh38/selfchain-LCR.hg38.bed \
      --out-dir "$GIAB/hg002_run/eval_v4_$SUB" || true
done
echo "DONE $(date)"
