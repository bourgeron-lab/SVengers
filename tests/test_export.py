import csv
import json

import numpy as np

from vedette.export import SCHEMA_VERSION, impact_fields
from vedette.report import _FamilySpill, _write_all_family_events


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


def test_family_export_rows(tmp_path):
    """genetics-viz issues #5/#7/#10: the family export carries every row of
    the events page, typed and coded columns, 1-based inclusive bounds, the
    per-child inheritance with the de novo DNQ, and the members' QC."""
    samples = ["kid", "dad", "mum"]
    sites = [
        {"id": "s1", "contig": "chr1", "start": 999, "end": 2000, "svtype": "DEL", "size": 1001,
         "filter": "PASS", "rescue": "", "idx": 0, "sqmax": 12.0, "carrier_frac": 0.001,
         "ref_carrier_frac": None},
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
    rows = {r["row_id"]: r for r in _rows(fam / "F1.vedette.tsv")}
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
    qc = {q["sample"]: q for q in _rows(fam / "F1.vedette.qc.tsv")}
    assert qc["kid"]["flags"] == "carrier_load;denovo_outlier;wavy"
    assert qc["kid"]["pass_sites_carried"] == "19626" and qc["dad"]["flags"] == ""
    meta = json.loads((fam / "F1.vedette.json").read_text())
    assert meta["schema_version"] == SCHEMA_VERSION and meta["members"] == samples
    assert meta["cohort"] == "test"
