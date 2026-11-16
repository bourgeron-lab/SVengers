import csv
import json

import numpy as np

from svengers.export import SCHEMA_VERSION, impact_fields, read_family_meta
from svengers.report import _FamilySpill, _write_all_family_events


def test_impact_fields():
    assert impact_fields("+0.30 — faible · DEL exonique SAMD11 (NDD candidat, s_het 0.01)") == (
        "0.30", "low", "DEL exonique SAMD11 (NDD candidat, s_het 0.01)")
    assert impact_fields("+1.90 — élevé · DEL complète SHANK3 (NDD, dominant) · 12 gènes codants")[:2] == (
        "1.90", "high")
    assert impact_fields("-0.60 — aucun gène codant") == ("-0.60", "none", "")
    assert impact_fields("") == ("", "", "")


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _own_call_vcf(path, sample, records):
    """A per-sample filtered VCF: (contig, start0, end, svtype, filter, source, imprecise, fmt)."""
    import pysam

    path.parent.mkdir(parents=True, exist_ok=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##svengers_coords=padded")
    hdr.contigs.add("chr1", length=1_000_000)
    for k, d in (("DEL", "Deletion"), ("DUP", "Duplication")):
        hdr.add_line(f'##ALT=<ID={k},Description="{d}">')
    hdr.add_line('##FILTER=<ID=LowSQ,Description="low SQ">')
    for k, n, t in (("SVTYPE", "1", "String"), ("END", "1", "Integer"), ("SOURCE", "1", "String"),
                    ("CIPOS", "2", "Integer"), ("CIEND", "2", "Integer"),
                    ("IMPRECISE", "0", "Flag"), ("MODEL_VER", "1", "String")):
        hdr.info.add(k, n, t, k)
    for k, t in (("GT", "String"), ("GQ", "Integer"), ("SQ", "Float"), ("PE", "Integer"),
                 ("SR", "Integer"), ("CG", "Integer"), ("RD", "Float"), ("CNQ", "Float")):
        hdr.formats.add(k, "1", t, k)
    hdr.add_sample(sample)
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for contig, start0, end, svtype, filt, source, imprecise, fmt in records:
            rec = vf.new_record(contig=contig, start=start0 - 1, stop=end,
                                alleles=("N", f"<{svtype}>"), id=f"svengers_{svtype}_{start0}")
            rec.info["SVTYPE"], rec.info["SOURCE"], rec.info["MODEL_VER"] = svtype, source, "v5n"
            rec.info["CIPOS"], rec.info["CIEND"] = (-5, 5), (-10, 10)
            if imprecise:
                rec.info["IMPRECISE"] = True
            rec.filter.add(filt)
            for k, v in fmt.items():
                rec.samples[sample][k] = v
            vf.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_family_export_rows(tmp_path):
    """genetics-viz issues #5/#7/#10, #8, #9, #6: the family export carries
    every row of the events page, typed and coded columns, 1-based inclusive
    bounds, the per-child inheritance with the de novo DNQ, the members' QC,
    each member's own call (SOURCE, precision, SQ/PE/SR...), the cohort
    site's key and evidence, and the same rows in long Parquet form."""
    samples = ["kid", "dad", "mum"]
    _own_call_vcf(tmp_path / "samples" / "kid" / "kid.filtered.vcf.gz", "kid", [
        ("chr1", 999, 2000, "DEL", "PASS", "BOTH", False,
         {"GT": (0, 1), "GQ": 60, "SQ": 12.3, "PE": 5, "SR": 3, "CG": 0, "RD": 0.52, "CNQ": 0.98}),
        ("chr1", 9999, 12000, "DUP", "LowSQ", "DEPTH", True,
         {"GT": (0, 1), "GQ": 30, "SQ": 2.1, "PE": 0, "SR": 0, "CG": 0, "RD": 1.45, "CNQ": 0.7}),
    ])
    sites = [
        {"id": "s1", "contig": "chr1", "start": 999, "end": 2000, "svtype": "DEL", "size": 1001,
         "filter": "PASS", "rescue": "", "idx": 0, "sqmax": 12.0, "carrier_frac": 0.001,
         "ref_carrier_frac": None, "af": 0.0005, "hwe": 1.0, "nsamples": 1,
         "precise_frac": 1.0, "sqmean": 12.0, "cn_disc": 0.0, "ref_af": 0.0},
        {"id": "s2", "contig": "chr1", "start": 9999, "end": 12000, "svtype": "DUP", "size": 2001,
         "filter": "PASS", "rescue": "", "idx": 1, "sqmax": None, "carrier_frac": None,
         "ref_carrier_frac": None},
    ]
    nalt = np.array([[1, 0, 0], [1, 1, 0]], dtype=np.int8)    # s1 de novo, s2 paternal
    gq = np.array([[60, 99, 99], [99, 60, 99]], dtype=np.int16)
    cn = np.array([[1, 2, 2], [3, 3, 2]], dtype=np.int8)
    inh = _FamilySpill(tmp_path / "inh")
    inh.add("F1", "s1\tkid\tdenovo\t45")
    inh.add("F1", "s2\tkid\tpaternal\t")

    def impact(contig, start, end, svtype="DEL"):
        return 0.9, "+0.90 — élevé · DEL exonique GENE1 (NDD, dominant)"

    def gene_lists(contig, start, end):
        return ([("GENE1", "ENSG1")], [("GENE2", "ENSG2")]) if start < 2000 else ([], [])

    info = {"kid": {"pass_sites_carried": 19626, "carrier_excess": 5.47, "denovo_candidates": 460,
                    "flags": ["carrier_load", "denovo_outlier"]}}
    n = _write_all_family_events(tmp_path, {"F1": samples}, samples, sites, nalt, gq, cn,
                                 lambda *a: "GENE1", impact, inh,
                                 sample_qc={"kid": ["charge 5.5×", "de novo 460"]},
                                 wave_scores={"kid": 0.2}, gene_lists_fn=gene_lists,
                                 sample_info=info, export_meta={"cohort": "test"})
    assert n == 1
    fam = tmp_path / "families" / "F1"
    rows = {r["row_id"]: r for r in _rows(fam / "F1.svengers.tsv")}
    r = rows["s1"]
    assert (r["kind"], r["site_ids"], r["start"], r["end"], r["vcf_pos"], r["size"]) == (
        "site", "s1", "1000", "2000", "999", "1001")
    assert (r["impact_score"], r["impact_level"]) == ("0.90", "high")
    assert r["impact_reason"] == "DEL exonique GENE1 (NDD, dominant)"
    assert (r["exonic_genes"], r["exonic_gene_ids"], r["genic_genes"]) == ("GENE1", "ENSG1", "GENE2")
    assert (r["GT_kid"], r["GQ_kid"], r["CN_kid"], r["INH_kid"], r["DNQ_kid"]) == (
        "0/1", "60", "1", "denovo", "45")
    assert (r["GT_dad"], r["INH_dad"]) == ("0/0", "")
    assert r["carriers"] == "kid" and r["sqmax"] == "12.0"
    assert "charge 5.5×" in r["qc_alerts"]
    assert (rows["s2"]["INH_kid"], rows["s2"]["DNQ_kid"]) == ("paternal", "")
    assert r["inheritance"] == "kid:denovo" and r["genes"] == "GENE1"
    qc = {q["sample"]: q for q in _rows(fam / "F1.svengers.qc.tsv")}
    assert qc["kid"]["flags"] == "carrier_load;denovo_outlier;wavy"
    assert qc["kid"]["pass_sites_carried"] == "19626" and qc["dad"]["flags"] == ""
    meta = json.loads((fam / "F1.svengers.json").read_text())
    assert meta["schema_version"] == SCHEMA_VERSION and meta["members"] == samples
    assert meta["cohort"] == "test"
    assert read_family_meta(fam, "F1")["svengers_version"] == meta["svengers_version"]

    # #9: the cohort site's own bounds; #8: site evidence and own calls
    assert r["site_key"] == "chr1:1000-2000:del" and rows["s2"]["site_key"] == "chr1:10000-12000:dup"
    assert (r["af"], r["nsamples"], r["precise_frac"], r["hwe_p"]) == ("0.00050", "1", "1.000", "1")
    assert (r["OWN_kid"], r["SQ_kid"], r["SRC_kid"], r["PRECISE_kid"]) == ("pass", "12.30", "BOTH", "1")
    assert (r["OWN_dad"], r["SQ_dad"]) == ("none", "")
    s2 = rows["s2"]
    assert (s2["OWN_kid"], s2["SRC_kid"], s2["PRECISE_kid"]) == ("filtered", "DEPTH", "0")

    # #6: long Parquet, one line per row x member, typed, versions in the metadata
    import pyarrow.parquet as pq

    assert meta["parquet"] == "F1.svengers.parquet"
    t = pq.read_table(fam / "F1.svengers.parquet")
    assert t.num_rows == 2 * len(samples)
    assert str(t.schema.field("start").type) == "int64" and str(t.schema.field("sq").type) == "double"
    import pyarrow as pa

    assert pa.types.is_list(t.schema.field("exonic_genes").type)
    md = {k.decode(): v.decode() for k, v in t.schema.metadata.items()}
    assert md["schema_version"] == str(SCHEMA_VERSION) and md["cohort"] == "test"
    assert json.loads(md["svengers"])["members"] == samples
    d = t.to_pylist()
    kid = next(x for x in d if x["row_id"] == "s1" and x["sample"] == "kid")
    assert (kid["carrier"], kid["gt"], kid["cn"], kid["inheritance"], kid["dnq"]) == (
        True, "0/1", 1, "denovo", 45)
    assert (kid["own_call_id"], kid["own_filter"], kid["source"], kid["precise"]) == (
        "svengers_DEL_999", "PASS", "BOTH", True)
    assert (kid["own_start"], kid["own_end"], kid["cipos"], kid["pe"], kid["sr"]) == (
        1000, 2000, [-5, 5], 5, 3)
    assert kid["sq"] == 12.3 or abs(kid["sq"] - 12.3) < 1e-4
    assert kid["own_model_version"] == "v5n" and kid["exonic_genes"] == ["GENE1"]
    dad = next(x for x in d if x["row_id"] == "s1" and x["sample"] == "dad")
    assert (dad["carrier"], dad["member_own_call"], dad["own_filter"], dad["sq"]) == (
        False, "none", None, None)
