#!/usr/bin/env bash
# Build the GRCh38 resource directory that configs/grch38.toml points to, from public
# sources, and check every file against the checksum of the copy the paper results were
# produced with (model v12b was trained on these exact files: a different version of a
# track changes the filter's features and its scores).
#
#   workflows/build_grch38_resources.sh OUT_DIR CLINSV_REFDATA
#
# OUT_DIR         directory to fill (created; existing files are re-checked, not re-made)
# CLINSV_REFDATA  ClinSV GRCh38 reference data v1.0, either the extracted `refdata-b38/`
#                 directory or the `refdata-b38_v1.0.tar` archive (38.95 GB), downloaded by
#                 you from https://clinsv.s3.ccia.org.au/clinsv_b38/refdata-b38_v1.0.tar
#                 under the ClinSV licence (research and education use; no redistribution,
#                 which is why this script cannot ship or fetch it for you). Three tracks
#                 come from it: popCovStdev.bw, the segmental duplications and the CNVnator
#                 exclusion list.
#
# Downloads (about 1.6 GB): the GIAB GRCh38 analysis-set FASTA (791 MB), the UCSC hg38
# Self Chain table (776 MB) and LCR-hs38 (12 MB). Needs: curl, gzip, awk, sort, tar,
# bedtools, md5sum, and Python with pysam (the SVengers environment) for the FASTA index.
set -euo pipefail

OUT=${1:?usage: build_grch38_resources.sh OUT_DIR CLINSV_REFDATA}
CLINSV=${2:?usage: build_grch38_resources.sh OUT_DIR CLINSV_REFDATA}
PYTHON=${PYTHON:-python3}
export LC_ALL=C
mkdir -p "$OUT"
cd "$OUT"

FASTA=GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
GIAB_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/$FASTA.gz
SELFCHAIN_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chainSelf.txt.gz
LCR_URL=https://raw.githubusercontent.com/lh3/varcmp/master/scripts/LCR-hs38.bed.gz

# checksums of the files behind the paper (SVengers v2.1.0, model v12b)
cat > MD5SUMS.expected <<EOF
99be410b684e3fdb04ba67d98adba249  $FASTA
a5bf6abc1f5efecc0a510f3120d08a3b  $FASTA.fai
b4406771c8fe431d114430280b39bb60  selfchain-LCR.hg38.bed
cd033c2a4030fc94d62f3123b87bd796  exclude.cnvnator_100bp.GRCh38.20170403.bed
b848173719156b7a4daa71048e0100b0  GRCh38GenomicSuperDup.bed.gz
33e33b4d925da7a8dd1b55e25aa76504  popCovStdev.bw
EOF

log() { echo "[$(date +%H:%M:%S)] $*"; }

# 1. reference FASTA (GIAB GRCh38 v3 analysis set, false duplications masked) + index
if [ ! -s "$FASTA" ]; then
  log "downloading the GIAB GRCh38 FASTA"
  curl -fsSL "$GIAB_URL" | gzip -dc > "$FASTA.part" && mv "$FASTA.part" "$FASTA"
fi
[ -s "$FASTA.fai" ] || "$PYTHON" -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$FASTA"

# 2. low-complexity + self-chain mask, the slivar recipe for hg38
#    (https://github.com/brentp/slivar/wiki/finding-deletions-in-parent-child-duos):
#    UCSC Self Chain with normScore > 90, target side and query side; the query
#    coordinates are taken as stored (no minus-strand conversion) - kept as is, the
#    model learned on this file - plus LCR-hs38, contigs without "chr", merged.
if [ ! -s selfchain-LCR.hg38.bed ]; then
  log "building selfchain-LCR.hg38.bed (streams the 776 MB UCSC Self Chain table)"
  ( curl -fsSL "$SELFCHAIN_URL" | gzip -dc \
      | awk 'BEGIN{FS=OFS="\t"} $13 > 90 { print $3,$5,$6; print $7,$10,$11 }';
    curl -fsSL "$LCR_URL" | gzip -dc ) \
    | sed -e 's/^chr//' | sort -k1,1 -k2,2n -k3,3n | bedtools merge > selfchain-LCR.hg38.bed.part
  mv selfchain-LCR.hg38.bed.part selfchain-LCR.hg38.bed
fi

# 3. the three ClinSV tracks
CLINSV_FILES=(exclude.cnvnator_100bp.GRCh38.20170403.bed tracks/GRCh38GenomicSuperDup.bed.gz
              tracks/GRCh38GenomicSuperDup.bed.gz.tbi tracks/popCovStdev.bw)
if [ -d "$CLINSV" ]; then  # accept .../clinsv, .../clinsv/refdata-b38 or .../refdata-b38
  for d in "$CLINSV/refdata-b38" "$CLINSV/clinsv/refdata-b38" "$CLINSV"; do
    [ -e "$d/tracks/popCovStdev.bw" ] && { CLINSV=$d; break; }
  done
fi
for f in "${CLINSV_FILES[@]}"; do
  [ -s "$(basename "$f")" ] && continue
  if [ -d "$CLINSV" ]; then
    cp "$CLINSV/$f" .
  else
    log "extracting $f from $CLINSV (reads the whole archive)"
    tar -xOf "$CLINSV" "clinsv/refdata-b38/$f" > "$(basename "$f")"
  fi
done

# 4. check
log "checking against the paper's files"
if md5sum -c MD5SUMS.expected; then
  log "OK: $OUT matches the resources of the paper; point configs/grch38.toml at it"
else
  log "MISMATCH: a source has changed since the paper; results will differ from the published ones"
  exit 1
fi
