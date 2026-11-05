"""genetics-viz validations as labels (issue #11)."""

import csv

from vedette.filter.schema import FEATURES
from vedette.train.validations import (
    breakpoint_offsets,
    label_features,
    read_genetics_viz,
    write_labels,
)

HEADER = ["FID", "Variant", "Sample", "User", "Inheritance", "Validation", "CuratedStart",
          "CuratedEnd", "Comment", "Ignore", "Timestamp", "SiteId", "CallerVersion", "ModelVersion"]


def _validations(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        for r in rows:
            w.writerow(r)


def test_genetics_viz_validations_become_labels(tmp_path):
    """Append-only file: the last row of a (FID, Variant, Sample) not marked
    Ignore wins; present/absent are 1/0, uncertain and different are kept but
    are no training label; Variant and curated bounds are 1-based inclusive."""
    p = tmp_path / "vedette.tsv"
    v = "vedette-cohort-1.7.0"
    _validations(p, [
        ["F1", "chr1:10001-20000:del", "kid", "u1", "de novo", "uncertain", "", "", "", "0",
         "2026-11-01T10:00:00", "vedette_site_1", v, "v5n"],
        ["F1", "chr1:10001-20000:del", "kid", "u1", "de novo", "present", "10051", "19950", "SR", "0",
         "2026-11-02T10:00:00", "vedette_site_1", v, "v5n"],
        ["F1", "chr1:50001-60000:dup", "kid", "u1", "", "absent", "", "", "", "0",
         "2026-11-02T11:00:00", "vedette_site_2", v, "v5n"],
        ["F1", "chr1:50001-60000:dup", "kid", "u2", "", "present", "", "", "mistake", "1",
         "2026-11-03T09:00:00", "vedette_site_2", v, "v5n"],          # ignored
        ["F1", "chr2:1001-9000:dup", "kid", "u1", "", "different", "", "", "DEL", "0",
         "2026-11-02T12:00:00", "vedette_site_3", v, "v5n"],
        ["F1", "chr3:5001-8000:del", "dad", "u1", "", "present", "", "", "", "0",
         "2026-11-02T12:00:00", "vedette_site_4", v, "v11"],
        ["F1", "nonsense", "kid", "u1", "", "present", "", "", "", "0", "", "", v, "v5n"],
    ])
    vals, counts = read_genetics_viz(p)
    assert counts == {"rows": 7, "ignored": 1, "superseded": 1, "unreadable": 1}
    by = {(x.contig, x.start): x for x in vals}
    assert len(vals) == 4
    first = by[("chr1", 10_000)]
    assert (first.end, first.svtype, first.validation, first.label) == (20_000, "DEL", "present", 1)
    assert (first.curated_start, first.curated_end) == (10_050, 19_950)
    assert by[("chr1", 50_000)].label == 0
    assert by[("chr2", 1_000)].label is None
    assert breakpoint_offsets(vals) == {"start": [50], "end": [-50]}

    out = tmp_path / "labels.tsv"
    write_labels(vals, out)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["label"] for r in rows if r["sample"] == "kid"] == ["1", "0", ""]
    kid = next(r for r in rows if r["sample"] == "kid" and r["start"] == "10000")
    assert kid["start_offset"] == "50" and kid["end_offset"] == "-50"

    # training: one sample's feature rows, matched by type and reciprocal overlap
    feats = tmp_path / "kid.features.tsv"
    cols = ["call_id", "contig", "start", "end", "svtype", *FEATURES]
    with open(feats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for cid, contig, s, e, t in [("c1", "chr1", 10_100, 19_900, "DEL"),   # present
                                     ("c2", "chr1", 50_000, 60_000, "DUP"),   # absent
                                     ("c3", "chr1", 50_000, 60_000, "DEL"),   # other type
                                     ("c4", "chr2", 1_000, 9_000, "DUP"),     # different
                                     ("c5", "chr3", 5_000, 8_000, "DEL")]:    # dad's
            w.writerow([cid, contig, s, e, t, *([0] * len(FEATURES))])
    labeled = tmp_path / "kid.labeled.tsv"
    c = label_features(feats, vals, "kid", labeled)
    assert c == {"tp": 1, "fp": 1, "unlabeled": 3, "conflict": 0}
    got = [(r["call_id"], r["label"]) for r in csv.DictReader(open(labeled), delimiter="\t")]
    assert got == [("c1", "1"), ("c2", "0")]


def test_validations_without_version_columns(tmp_path):
    """svs.tsv-like files stop at Timestamp: the version columns are optional."""
    p = tmp_path / "vedette.tsv"
    with open(p, "w") as fh:
        fh.write("\t".join(HEADER[:11]) + "\n")
        fh.write("F1\tchr1:1,001-5,000:DEL\tkid\tu\t\tpresent\t\t\t\t0\t2026-11-01T10:00:00\n")
    vals, _ = read_genetics_viz(p)
    assert [(x.start, x.end, x.svtype, x.model_version) for x in vals] == [(1000, 5000, "DEL", "")]
