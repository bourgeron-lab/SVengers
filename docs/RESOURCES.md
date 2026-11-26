# GRCh38 resources

SVengers needs a reference FASTA and four companion tracks. The filter models
(`models/v15`, `models/v12b`) were trained on features computed from these exact files: a
different version of a track changes the features and the scores, so results
match the published ones only with the files below. They are not in the
repository; `workflows/build_grch38_resources.sh` rebuilds the directory from
the public sources and checks each file against the checksum of the copy the
paper used. `configs/grch38.toml` then points at that directory.

```bash
workflows/build_grch38_resources.sh /path/to/svengers-grch38 /path/to/clinsv/refdata-b38
```

The second argument is the ClinSV GRCh38 reference data v1.0 (extracted
directory or the `refdata-b38_v1.0.tar` archive), which you download yourself:
its licence does not allow redistribution (see below). The script downloads
about 1.6 GB and needs curl, gzip, awk, sort, tar, bedtools, md5sum and Python
with pysam.

## Files

| config key | file | used by | source | md5 |
|---|---|---|---|---|
| `fasta` | `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta` (+ `.fai`, md5 `a5bf6abc1f5efecc0a510f3120d08a3b`) | every stage | GIAB GRCh38 v3 analysis set with the GRC false duplications masked, [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/) (`.fasta.gz`, 791 MB) | `99be410b684e3fdb04ba67d98adba249` |
| `lcr_bed` | `selfchain-LCR.hg38.bed` | filter (`lcr_ovl`), report | built, see below | `b4406771c8fe431d114430280b39bb60` |
| `exclude_bed` | `exclude.cnvnator_100bp.GRCh38.20170403.bed` | call (bin mask), filter (`exclude_ovl`) | ClinSV reference data; originally the speedseq annotation of the same name (same size, not checked byte for byte) | `cd033c2a4030fc94d62f3123b87bd796` |
| `segdup_bed` | `GRCh38GenomicSuperDup.bed.gz` (+ `.tbi`) | filter (`segdup_ovl`), report | ClinSV reference data (UCSC genomicSuperDups, IGV-formatted) | `b848173719156b7a4daa71048e0100b0` |
| `pop_cov_stdev_bw` | `popCovStdev.bw` | filter (`popcov_stdev`) | ClinSV reference data (`tracks/popCovStdev.bw`) | `33e33b4d925da7a8dd1b55e25aa76504` |

GC content and N fraction per bin are computed from the FASTA on first use and
cached in the `--resources` directory of `svengers call`.

### `selfchain-LCR.hg38.bed`

The recipe of the [slivar wiki](https://github.com/brentp/slivar/wiki/finding-deletions-in-parent-child-duos)
for hg38, rebuilt byte for byte on 2026-11-17:

```bash
( curl -fsSL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chainSelf.txt.gz | gzip -dc \
    | awk 'BEGIN{FS=OFS="\t"} $13 > 90 { print $3,$5,$6; print $7,$10,$11 }';
  curl -fsSL https://raw.githubusercontent.com/lh3/varcmp/master/scripts/LCR-hs38.bed.gz | gzip -dc ) \
  | sed -e 's/^chr//' | sort -k1,1 -k2,2n -k3,3n | bedtools merge > selfchain-LCR.hg38.bed
```

UCSC hg38 Self Chain chains with `normScore > 90` (table of 2020-06-21), both
target and query sides, plus Heng Li's low-complexity regions (LCR-hs38).
Two quirks are kept on purpose because the model learned on them: the query
coordinates are used as stored, without the minus-strand conversion, so about
41 Mb (12%) of the mask sit at mirrored positions; and contig names have no
`chr` prefix (SVengers matches both conventions).

### ClinSV tracks and licence

`exclude_bed`, `segdup_bed` and `pop_cov_stdev_bw` are the copies shipped in
the [ClinSV](https://github.com/KCCG/ClinSV) GRCh38 reference data v1.0
(`https://clinsv.s3.ccia.org.au/clinsv_b38/refdata-b38_v1.0.tar`, 38.95 GB).
The ClinSV licence grants research and education use and does not allow
redistribution, so these files cannot ship with SVengers or in a resource
archive. The exclusion list and the segmental duplications have public
origins (speedseq, UCSC); `popCovStdev.bw` is specific to ClinSV.

## Effect of a missing track

Measured on HG002 (GIAB v5.0q deletions, model v12b retrained without HG002)
and on 186 genomes of the paper cohort with their 20 trios
(`runs/_repro/` on the GHFC cluster):

| configuration | HG002 DEL F1 (all / 1-10 kb) | PASS per genome | trio children |
|---|---|---|---|
| all tracks | 0.515 / 0.836 | — | — |
| without `pop_cov_stdev_bw` | 0.516 / 0.801 | DEL +12% (noise < 1 kb), DUP -27% | ~80 real 1-10 kb deletions lost per child |
| without `lcr_bed` | 0.503 / 0.833 | DEL -12%, DUP +26% (noise < 1 kb) | 180 real 1-10 kb deletions lost in 20 children |

The code runs without either track (the feature is then 0), but the scores
are no longer those of the paper. A substitute track (another LCR mask,
another population-depth track) is a different feature and needs a retrained
model; this is planned after the paper (`docs/ROADMAP.md`).

## Not distributed

- **Panel of normals** (`--pon`): built from your own cohort, with the same
  sequencing and processing as the genomes to call (README, *Panel of
  normals*). The paper used 300 founders of its cohort.
- **Filter model**: `models/v15` (default since v2.3.0) and `models/v12b` (the
  paper's) are in the repository; their training data are described in
  `docs/MODEL.md`.
- **Ensembl GFF** (gene column of the report) and the **GeneTrek** gene table
  (impact score): optional.
