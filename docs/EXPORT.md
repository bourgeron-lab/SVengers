# Family export format

`vedette report --family-events` writes, next to each `families/<FID>/events.html`:

| file | content |
|---|---|
| `<FID>.vedette.tsv` | one row per row of the events page |
| `<FID>.vedette.qc.tsv` | one row per family member: the per-sample QC behind the row alerts |
| `<FID>.vedette.json` | `schema_version`, `vedette_version`, `model_version`, cohort, member order, options |

**`schema_version`** (currently **1**, in the JSON) changes whenever a column is
renamed, removed or changes meaning. Adding a column does not change it: read
columns by name.

**Missing values** are empty fields, everywhere. Lists are `;`-separated.

## Coordinates

- `start`, `end`: first and last **affected** base, 1-based, both included
  (what IGV and UCSC show). `size = end - start + 1`.
- `vcf_pos = start - 1`: the VCF `POS`, the padding base before the event.

In the VCFs (from v1.4.0, header line `##vedette_coords=padded`), records follow
the VCF 4.3 rule for symbolic alleles: `POS` is the base before the event,
`END` its last base, `SVLEN = ±(END - POS)`. A deletion whose clipped reads
show the last retained base at 1,000 and the first base after the gap at
1,101 is `POS=1000`, `END=1100`, `SVLEN=-100`, and `start=1001`, `end=1100` here.

Vedette VCFs written before v1.4.0 (no `##vedette_coords` line) had
`POS = start` and `END = end + 1` in the terms above, one base to the right on
both sides. Vedette reads both. Breakpoints curated against the old files are
off by one at the start and at the end.

## `<FID>.vedette.tsv`

Rows are the rows of the events page. Cohort sites carried by the same members
are collapsed into one event (`site_ids` lists them), and their bounds are the
carriers' own calls when there are any. Sub-threshold, mosaic and chromosomal
rows are added. With `--events-pass-only` (GHFC production), non-PASS cohort
sites are left out of the page and of this file.

| column | type | values |
|---|---|---|
| `row_id` | string | first cohort site id (`vedette_site_N`), or `sub_<S>_<chr>_<pos>`, `mosaic_…`, `chrom_…` |
| `kind` | code | `site`, `subthreshold`, `mosaic`, `chromosomal` |
| `site_ids` | list | cohort site ids behind a `site` row; empty otherwise |
| `contig` | string | |
| `start`, `end` | int | see Coordinates |
| `vcf_pos` | int | `start - 1` |
| `svtype` | code | `DEL`, `DUP` |
| `size` | int | `end - start + 1` |
| `filter` | string | cohort FILTER (`PASS`, or `;`-joined names); `LowSQ` for sub-threshold rows |
| `rescue` | string | `RECURRENT`, `SOUS-SEUIL NDD`, `SOUS-SEUIL GRANDE DUP` (from v1.6.0), `MOSAÏQUE`, `CHROMOSOMIQUE` or empty |
| `n_sites` | int | cohort sites collapsed into the row (`site` rows) |
| `sqmax` | float | best per-sample SQ among the site's calls (sub-threshold rows: the call's SQ) |
| `cohort_carrier_frac` | float | fraction of the cohort carrying the site |
| `ref_carrier_frac` | float | the same in the reference cohort (`--known-sites`) |
| `lcr_ovl`, `segdup_ovl` | float | fraction of the row in LCR / segmental duplications |
| `support` | code | `strong`: a carrier's re-genotyped copy number agrees with the type; `weak`: none does |
| `own_call` | code | `pass`: a member's own call (before merging) passed the filter; `filtered`: an own call exists but was filtered; `none`: carriage from cohort re-genotyping only; `lowsq` for sub-threshold rows; empty for mosaic/chromosomal rows |
| `mosaic_cn` | float | fractional copy number of mosaic and chromosomal rows |
| `impact_score` | float | gene-content score (ACMG/ClinGen CNV sections 1-3: ClinGen dosage regions and genes, GeneTrek for the genes ClinGen has not scored) |
| `impact_level` | code | `high` (≥ 0.90), `moderate` (≥ 0.45), `low` (≥ 0), `none` (no coding gene), `benign` (−1.00: inside a ClinGen region curated "dosage sensitivity unlikely" for its type) |
| `impact_reason` | text | the rule and gene behind the score (French) |
| `impact` | text | display label (French), the three fields above in one string |
| `exonic_genes`, `exonic_gene_ids` | list | protein-coding genes with an exon in the row (symbols, Ensembl ids), position order, complete |
| `genic_genes`, `genic_gene_ids` | list | protein-coding genes overlapped without an exon (intronic) |
| `carriers` | list | members carrying the row (GT ≥ 0/1 and GQ ≥ `min_carrier_gq`) |
| `qc_alerts` | text | the page's alert column (French) |
| `genes` | text | display list (coding genes, truncated after six names); kept from the pre-v1.4.0 export, prefer `exonic_genes` / `genic_genes` |
| `inheritance` | text | `<child>:<code>; ...`, the `INH_<id>` columns in one field; kept from the pre-v1.4.0 export |

Then, for each member in the JSON's `members` order:

| column | type | values |
|---|---|---|
| `GT_<id>` | code | `0/0`, `0/1`, `1/1`, `./.`; empty on mosaic/chromosomal rows and for members without a call on sub-threshold rows |
| `GQ_<id>` | int | genotype quality |
| `CN_<id>` | int | re-genotyped copy number |
| `INH_<id>` | code | children of a trio only: `denovo`, `denovo_lowconf` (a trio member's GQ under the de novo threshold), `paternal`, `maternal`, `biparental`, `none`, `karyotype_<K>` (a de novo call on the sex chromosome the child's karyotype `<K>` explains); over the sites of a collapsed row, the label most of them carry |
| `DNQ_<id>` | int | de novo quality of `denovo` rows: min GQ of the trio |

## `<FID>.vedette.qc.tsv`

| column | type | values |
|---|---|---|
| `sample` | string | |
| `wave_score` | float | coverage waviness (≥ 0.128: `wavy` flag) |
| `pass_sites_carried` | int | PASS cohort sites the sample carries |
| `carrier_excess` | float | the same divided by the cohort median |
| `denovo_candidates` | int | PASS de novo candidates (children) |
| `karyotype` | string | `XX`, `XY`, `XXY`, ... or `?` |
| `flags` | list | `carrier_load` (≥ 5× the cohort median), `denovo_outlier` (de novo candidates above the cohort cutoff: left out of the report's tables), `karyotype` (atypical), `wavy`, `big_cnv` (≥ 3 PASS events ≥ 1 Mb), `mosaic` (mosaic segments), `subthreshold_masked` (too many sub-threshold rows: masked) |
| `flags_text` | text | the same as shown on the page (French) |
