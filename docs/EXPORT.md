# Family export format

`svengers report --family-events` writes, next to each `families/<FID>/events.html`
(before v2.0, when SVengers was Vedette: `<FID>.vedette.*`, with `vedette_version` and
the Parquet key `vedette`; `svengers.export.read_family_meta` reads both):

| file | content |
|---|---|
| `<FID>.svengers.tsv` | one row per row of the events page |
| `<FID>.svengers.qc.tsv` | one row per family member: the per-sample QC behind the row alerts |
| `<FID>.svengers.json` | `schema_version`, `svengers_version`, `model_version`, cohort, panel of normals (`pon`), member order, options |
| `<FID>.svengers.parquet` | from v1.9.0, when pyarrow is installed: the rows of the TSV in long form, one line per row × member, typed (see below) |

**`schema_version`** (currently **1**, in the JSON) changes whenever a column is
renamed, removed or changes meaning. Adding a column does not change it: read
columns by name.

**Missing values** are empty fields, everywhere. Lists are `;`-separated.

## Coordinates

- `start`, `end`: first and last **affected** base, 1-based, both included
  (what IGV and UCSC show). `size = end - start + 1`.
- `vcf_pos = start - 1`: the VCF `POS`, the padding base before the event.

In the VCFs (from v1.4.0, header line `##svengers_coords=padded`, `##vedette_coords=padded`
before v2.0), records follow
the VCF 4.3 rule for symbolic alleles: `POS` is the base before the event,
`END` its last base, `SVLEN = ±(END - POS)`. A deletion whose clipped reads
show the last retained base at 1,000 and the first base after the gap at
1,101 is `POS=1000`, `END=1100`, `SVLEN=-100`, and `start=1001`, `end=1100` here.

Vedette VCFs written before v1.4.0 (no `##vedette_coords` line) had
`POS = start` and `END = end + 1` in the terms above, one base to the right on
both sides. SVengers reads both. Breakpoints curated against the old files are
off by one at the start and at the end.

## `<FID>.svengers.tsv`

Rows are the rows of the events page. Cohort sites carried by the same members
are collapsed into one event (`site_ids` lists them), and their bounds are the
carriers' own calls when there are any. Sub-threshold, mosaic and chromosomal
rows are added. With `--events-pass-only` (GHFC production) the page holds
PASS sites only; this file also holds, from v1.7.0, the non-PASS sites carried
by at most 20 % of the cohort (`--export-nonpass-max-carrier-frac`, written
as `export_nonpass_max_carrier_frac` in the `.json`), as `site` rows whose
`filter` is not `PASS`: a reader shows them behind a toggle. The non-PASS
sites carried more widely (dataset_paper: 167,560 of 190,367, each carried by
over half the cohort) are reference artefacts and are left out. Non-PASS
sites are collapsed among themselves, never with PASS sites, and one that a
sub-threshold, mosaic or chromosomal row already covers is not repeated.

| column | type | values |
|---|---|---|
| `row_id` | string | first cohort site id (`svengers_site_N`, `vedette_site_N` before v2.0), or `sub_<S>_<chr>_<pos>`, `mosaic_…`, `chrom_…` |
| `kind` | code | `site`, `subthreshold`, `mosaic`, `chromosomal` |
| `site_ids` | list | cohort site ids behind a `site` row; empty otherwise |
| `site_key` | string | from v1.9.0: `chr:start-end:del\|dup` of the cohort site `row_id`, its own bounds (1-based inclusive), not the row's; see Site keys |
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
| `af`, `nsamples`, `precise_frac`, `sqmean`, `hwe_p`, `cn_disc`, `ref_af` | float / int | from v1.9.0, the cohort site `row_id`'s INFO: allele frequency, distinct samples clustered at the merge, fraction of them with precise breakpoints, mean of their best SQ, founders' Hardy-Weinberg p, fraction of carriers whose copy number contradicts the type, frequency in the reference cohort; empty on other rows and when the cohort VCF lacks the field |
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
| `OWN_<id>` | code | from v1.9.0: the member's own call (its filtered VCF, before merging) behind the row, same type and reciprocal overlap ≥ 0.5, a PASS call first: `pass`, `filtered`, `none`; empty on mosaic/chromosomal rows |
| `SQ_<id>` | float | that call's SQ |
| `SRC_<id>` | code | that call's evidence source: `DEPTH`, `JUNCTION`, `BOTH`, `INSERT` |
| `PRECISE_<id>` | code | `1`: base-pair breakpoints, `0`: imprecise |

The other fields of the own call (id, FILTER, bounds, CIPOS/CIEND, GQ, PE, SR,
CG, RD, CNQ, model version) are in the Parquet file.

## `<FID>.svengers.qc.tsv`

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

## `<FID>.svengers.parquet`

From v1.9.0, when pyarrow is installed (`pip install svengers[parquet]`; the
`.json` gives the file name in `parquet`, null when it was not written). The
rows of `<FID>.svengers.tsv` in long form: **one line per row × family
member**, non-carriers included, sorted by `contig`, `start`. Row-level
columns repeat on each member's line, with the TSV's names and meanings, typed
(int64, float64, string, bool; nulls are real nulls) and with lists as list
columns: `row_id`, `kind`, `site_ids` (list), `site_key`, `contig`, `start`,
`end`, `vcf_pos`, `svtype`, `size`, `filter`, `rescue`, `n_sites`, `sqmax`,
`cohort_carrier_frac`, `ref_carrier_frac`, `af`, `nsamples`, `precise_frac`,
`sqmean`, `hwe_p`, `cn_disc`, `ref_af`, `lcr_ovl`, `segdup_ovl`, `support`,
`own_call`, `mosaic_cn`, `impact_score`, `impact_level`, `impact_reason`,
`exonic_genes`, `exonic_gene_ids`, `genic_genes`, `genic_gene_ids`,
`carriers` (lists), `qc_alerts`. The display columns `impact`, `genes`,
`inheritance` are left out.

Member columns:

| column | type | values |
|---|---|---|
| `sample` | string | member id |
| `carrier` | bool | in `carriers` |
| `gt`, `gq`, `cn` | string, int, int | `GT_`, `GQ_`, `CN_` of the TSV |
| `inheritance`, `dnq` | string, int | `INH_`, `DNQ_` |
| `member_own_call` | string | `OWN_`: `pass`, `filtered`, `none` |
| `own_call_id`, `own_filter` | string | the own call's VCF id and FILTER |
| `own_start`, `own_end` | int | its bounds, 1-based inclusive |
| `source`, `precise` | string, bool | `SRC_`, `PRECISE_` |
| `cipos`, `ciend` | list of int | its CIPOS / CIEND |
| `sq`, `own_gq` | float, int | its SQ and GQ |
| `pe`, `sr`, `cg` | int | discordant pairs, split/clipped reads, reads with the event as a CIGAR deletion |
| `rd`, `cnq` | float | mean copy ratio over the call, copy-number posterior |
| `own_model_version` | string | filter model that scored the call (`MODEL_VER`) |

File metadata (Parquet key/value): `svengers` (the `.json`'s content, JSON),
and `schema_version`, `svengers_version`, `model_version`, `cohort`, `pon`,
`family_id`, `generated` as strings.

## Site keys and re-merges

A cohort site id (`svengers_site_N`) is renumbered at every merge, and the
merge moves a site's bounds when it sees other breakpoints. `site_key` is
stable while the site keeps its bounds; across a re-merge, `svengers
site-map OLD.vcf.gz NEW.vcf.gz --out map.tsv` gives the correspondence:

| column | values |
|---|---|
| `old_site_id`, `old_site_key` | the previous run's site (empty for `new`) |
| `new_site_id`, `new_site_key` | the new run's site (empty for `lost`) |
| `relation` | `same` (same bounds), `moved` (mutual best match by reciprocal overlap ≥ 0.5, other bounds), `split` (the old site matches several new ones), `merged` (several old sites match the new one), `complex` (both), `lost`, `new` |
| `overlap` | reciprocal overlap of the pair |

A validation keyed on an old `site_key` moves to the `new_site_key`(s) of its
`same`, `moved` or `merged` line; a `split` gives several candidates; `lost`
has none.
