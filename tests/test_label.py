"""Integration test of the Truvari labeling path with tiny synthetic VCFs.

Exercises the riskiest interface: Truvari accepting Vedette's symbolic
<DEL>/<DUP> records and --dup-to-ins matching.
"""

import shutil
import sys
from pathlib import Path

import pysam
import pytest

from vedette.filter.schema import FEATURES
from vedette.train.label import label_features, run_truvari_bench

truvari_missing = (
    shutil.which("truvari") is None
    and not Path(sys.executable).with_name("truvari").exists()
)


def _header(with_symbolic=False):
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.info.add("SVLEN", "1", "Integer", "l")
    h.info.add("END", "1", "Integer", "e")
    h.formats.add("GT", "1", "String", "g")
    if with_symbolic:
        h.add_line('##ALT=<ID=DEL,Description="Deletion">')
        h.add_line('##ALT=<ID=DUP,Description="Duplication">')
    h.add_sample("S")
    return h


def _write(path, header, records):
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, stop, alleles, svtype, svlen, rid in records:
            rec = out.new_record(
                contig=contig, start=start, stop=stop, alleles=alleles, id=rid
            )
            rec.info["SVTYPE"] = svtype
            rec.info["SVLEN"] = svlen
            rec.samples["S"]["GT"] = (0, 1)
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


@pytest.mark.skipif(truvari_missing, reason="truvari not installed")
def test_truvari_labeling_roundtrip(tmp_path):
    # truth: explicit-sequence DEL (500 bp) and an INS (600 bp) as DUP proxy
    truth = tmp_path / "truth.vcf.gz"
    _write(truth, _header(), [
        ("chr1", 100_000, 100_501, ("A" + "T" * 500, "A"), "DEL", -500, "t1"),
        ("chr1", 800_000, 800_001, ("A", "A" + "C" * 600), "INS", 600, "t2"),
    ])
    # calls: symbolic; c1 matches t1, c_dup matches t2 via --dup-to-ins,
    # c_fp matches nothing
    calls = tmp_path / "calls.vcf.gz"
    _write(calls, _header(with_symbolic=True), [
        ("chr1", 100_010, 100_505, ("N", "<DEL>"), "DEL", -495, "c1"),
        ("chr1", 500_000, 501_000, ("N", "<DEL>"), "DEL", -1000, "c_fp"),
        ("chr1", 799_990, 800_590, ("N", "<DUP>"), "DUP", 600, "c_dup"),
    ])
    bed = tmp_path / "conf.bed"
    bed.write_text("chr1\t0\t1000000\n")

    tdir = run_truvari_bench(calls, truth, bed, None, tmp_path / "bench")

    feats = tmp_path / "features.tsv"
    with open(feats, "w") as fh:
        cols = ["call_id", "contig", "start", "end", "svtype", *FEATURES]
        fh.write("\t".join(cols) + "\n")
        for rid, svtype in (("c1", "DEL"), ("c_fp", "DEL"), ("c_dup", "DUP")):
            fh.write("\t".join([rid, "chr1", "0", "1", svtype] + ["0"] * len(FEATURES)) + "\n")

    out = tmp_path / "labeled.tsv"
    counts = label_features(feats, tdir, out)
    labels = {
        line.split("\t")[0]: line.strip().split("\t")[-1]
        for line in out.read_text().splitlines()[1:]
    }
    assert labels["c1"] == "1"
    assert labels["c_fp"] == "0"
    assert labels["c_dup"] == "1", "--dup-to-ins should match the DUP to the truth INS"
    assert counts["tp"] == 2 and counts["fp"] == 1


@pytest.mark.skipif(truvari_missing, reason="truvari not installed")
def test_parent_mode_labeling(tmp_path):
    """Parent mode: truth match -> TP; non-match -> FP only on recurrent-FP loci."""
    truth = tmp_path / "truth.vcf.gz"
    _write(truth, _header(), [
        ("chr1", 100_000, 100_501, ("A" + "T" * 500, "A"), "DEL", -500, "t1"),
    ])
    calls = tmp_path / "calls.vcf.gz"
    _write(calls, _header(with_symbolic=True), [
        ("chr1", 100_010, 100_505, ("N", "<DEL>"), "DEL", -495, "p_tp"),      # transmis
        ("chr1", 500_000, 501_000, ("N", "<DEL>"), "DEL", -1000, "p_artifact"),  # FP récurrent
        ("chr1", 700_000, 700_400, ("N", "<DEL>"), "DEL", -400, "p_private"),    # privé -> drop
    ])
    # loci FP récurrents (les FP du proband)
    fp_vcf = tmp_path / "hg002_fp.vcf.gz"
    _write(fp_vcf, _header(with_symbolic=True), [
        ("chr1", 499_900, 501_100, ("N", "<DEL>"), "DEL", -1200, "fp1"),
    ])
    bed = tmp_path / "conf.bed"
    bed.write_text("chr1\t0\t1000000\n")

    tdir = run_truvari_bench(calls, truth, bed, None, tmp_path / "bench")
    feats = tmp_path / "features.tsv"
    with open(feats, "w") as fh:
        cols = ["call_id", "contig", "start", "end", "svtype", *FEATURES]
        fh.write("\t".join(cols) + "\n")
        for rid, start, end in (("p_tp", 100_010, 100_505),
                                ("p_artifact", 500_000, 501_000),
                                ("p_private", 700_000, 700_400)):
            fh.write("\t".join([rid, "chr1", str(start), str(end), "DEL"]
                               + ["0"] * len(FEATURES)) + "\n")

    out = tmp_path / "labeled.tsv"
    counts = label_features(
        feats, tdir, out, parent_mode=True, recurrent_fp_vcf=fp_vcf
    )
    labels = {
        line.split("\t")[0]: line.strip().split("\t")[-1]
        for line in out.read_text().splitlines()[1:]
    }
    assert labels == {"p_tp": "1", "p_artifact": "0"}
    assert counts == {"tp": 1, "fp": 1, "ambiguous": 1, "unlabeled": 0}


def test_label_by_transmission(tmp_path):
    import csv

    from vedette.train.label_trio import label_by_transmission

    def write_gt(path, rows):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(
                fh, fieldnames=["site_id", "contig", "start", "end", "svtype",
                                "sample", "gt", "gq"], delimiter="\t")
            w.writeheader()
            for site_id, c, s, e, t, gt, gq in rows:
                w.writerow({"site_id": site_id, "contig": c, "start": s, "end": e,
                            "svtype": t, "sample": "P", "gt": gt, "gq": gq})

    # father carries site A; both parents confident 0/0 at site B; site C low-GQ
    write_gt(tmp_path / "father.gt.tsv", [
        ("A", "chr1", 10_000, 20_000, "DEL", "0/1", 99),
        ("B", "chr1", 50_000, 60_000, "DEL", "0/0", 80),
        ("C", "chr1", 90_000, 95_000, "DUP", "0/0", 5),
    ])
    write_gt(tmp_path / "mother.gt.tsv", [
        ("A", "chr1", 10_000, 20_000, "DEL", "0/0", 90),
        ("B", "chr1", 50_000, 60_000, "DEL", "0/0", 75),
        ("C", "chr1", 90_000, 95_000, "DUP", "0/0", 99),
    ])
    with open(tmp_path / "child.features.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["call_id", "contig", "start", "end",
                                           "svtype", "dhffc"], delimiter="\t")
        w.writeheader()
        rows = [
            ("c1", "chr1", 10_100, 19_800, "DEL", 0.5),   # transmis (père) -> TP
            ("c2", "chr1", 50_200, 59_500, "DEL", 0.6),   # 0/0 confiant x2 -> FP
            ("c3", "chr1", 90_000, 95_000, "DUP", 1.5),   # père GQ 5 -> ambigu
            ("c4", "chr1", 500_000, 505_000, "DEL", 0.5), # aucun site -> ambigu
        ]
        for r in rows:
            w.writerow(dict(zip(w.fieldnames, r, strict=True)))

    counts = label_by_transmission(
        tmp_path / "child.features.tsv",
        [tmp_path / "father.gt.tsv", tmp_path / "mother.gt.tsv"],
        tmp_path / "child.labeled.tsv",
    )
    assert counts == {"tp": 1, "fp": 1, "ambiguous": 2}
    labeled = list(csv.DictReader(open(tmp_path / "child.labeled.tsv"), delimiter="\t"))
    assert {(r["call_id"], r["label"]) for r in labeled} == {("c1", "1"), ("c2", "0")}
