"""ped-annotate after v2.5 (PASS calls keep a non-reference genotype): a
parent's 0/0 is read at min(GQ, RQ), and the cohort-fraction allowance on
carriers outside the family holds for large sites only."""

import csv

import pysam

from svengers.cohort.mendel import run_ped_annotate
from svengers.ped import parse_ped

SAMPLES = ("C", "F", "M", "U1", "U2", "U3")


def _cohort(path, start, end, father_rq=None, unrelated=0):
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.formats.add("GT", "1", "String", "g")
    h.formats.add("GQ", "1", "Integer", "q")
    h.formats.add("RQ", "1", "Integer", "r")
    for s in SAMPLES:
        h.add_sample(s)
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        r = out.new_record(contig="chr1", start=start - 1, stop=end, alleles=("N", "<DEL>"), id="s")
        r.info["SVTYPE"] = "DEL"
        for smp in SAMPLES:
            carrier = smp == "C" or (smp.startswith("U") and int(smp[1]) <= unrelated)
            r.samples[smp]["GT"] = (0, 1) if carrier else (0, 0)
            r.samples[smp]["GQ"] = 50
            r.samples[smp]["RQ"] = father_rq if (smp == "F" and father_rq is not None) else 60
        out.write(r)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _own(path, start, end):
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.formats.add("GT", "1", "String", "g")
    h.add_sample("C")
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        r = out.new_record(contig="chr1", start=start - 1, stop=end, alleles=("N", "<DEL>"))
        r.info["SVTYPE"] = "DEL"
        r.samples["C"]["GT"] = (0, 1)
        out.write(r)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _run(tmp_path, start, end, **kw):
    cohort_kw = {k: kw.pop(k) for k in ("father_rq", "unrelated") if k in kw}
    _cohort(tmp_path / "cohort.vcf.gz", start, end, **cohort_kw)
    d = tmp_path / "samples" / "C"
    d.mkdir(parents=True, exist_ok=True)
    _own(d / "C.filtered.vcf.gz", start, end)
    (tmp_path / "ped.tsv").write_text(
        "FAM\tF\t0\t0\t1\t1\nFAM\tM\t0\t0\t2\t1\nFAM\tC\tF\tM\t1\t2\n"
        + "".join(f"O{i}\tU{i}\t0\t0\t1\t1\n" for i in (1, 2, 3)))
    run_ped_annotate(tmp_path / "cohort.vcf.gz", parse_ped(tmp_path / "ped.tsv"),
                     tmp_path / "m.tsv", None, calls_dir=tmp_path / "samples", **kw)
    [row] = list(csv.DictReader(open(tmp_path / "m.tsv"), delimiter="\t"))
    return row


def test_parent_reference_quality_caps_confidence(tmp_path):
    assert _run(tmp_path, 10_000, 10_200)["inheritance"] == "denovo"
    row = _run(tmp_path, 10_000, 10_200, father_rq=10)
    assert row["inheritance"] == "denovo_lowconf" and row["father_gq"] == "10"


def test_carrier_fraction_holds_for_large_sites_only(tmp_path):
    # 3 carriers in other families; floor 2, fraction 0.5 of 6 samples = 3
    kw = dict(unrelated=3, max_other_carriers=2, max_other_carrier_frac=0.5)
    assert _run(tmp_path, 10_000, 10_500, **kw)["inheritance"] == "denovo_weak"
    assert _run(tmp_path, 10_000, 30_000, **kw)["inheritance"] == "denovo"
    assert _run(tmp_path, 10_000, 10_500, max_other_frac_min_size=0, **kw)["inheritance"] == "denovo"
