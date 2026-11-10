"""ped-annotate --calls-dir: a de novo without a call of the child at the site is genotype-only."""

import csv

import pysam

from vedette.cohort.mendel import OwnCalls, run_ped_annotate
from vedette.ped import parse_ped

SAMPLES = ("C", "F", "M")


def _cohort(path, sites):
    """sites: [(id, svtype, start, end, child_gt)], parents 0/0, GQ 50 everywhere."""
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.formats.add("GT", "1", "String", "g")
    h.formats.add("GQ", "1", "Integer", "q")
    for s in SAMPLES:
        h.add_sample(s)
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        for sid, t, s, e, cgt in sites:
            # spec coordinates: POS = padding base, END = last base
            r = out.new_record(contig="chr1", start=s - 1, stop=e, alleles=("N", f"<{t}>"), id=sid)
            r.info["SVTYPE"] = t
            for smp in SAMPLES:
                r.samples[smp]["GT"] = cgt if smp == "C" else (0, 0)
                r.samples[smp]["GQ"] = 50
            out.write(r)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _calls(path, calls):
    """calls: [(svtype, start, end, filter, gt, sq)] of sample C."""
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.filters.add("LowSQ", None, None, "low")
    h.formats.add("GT", "1", "String", "g")
    h.formats.add("SQ", "1", "Float", "s")
    h.add_sample("C")
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        for t, s, e, flt, gt, sq in calls:
            r = out.new_record(contig="chr1", start=s - 1, stop=e, alleles=("N", f"<{t}>"))
            r.info["SVTYPE"] = t
            if flt != "PASS":
                r.filter.add(flt)
            r.samples["C"]["GT"] = gt
            r.samples["C"]["SQ"] = sq
            out.write(r)
    pysam.tabix_index(str(path), preset="vcf", force=True)


SITES = [
    ("own_pass", "DEL", 10_000, 15_000, (0, 1)),
    ("own_lowsq", "DUP", 30_000, 32_000, (0, 1)),
    ("wrong_type", "DEL", 50_000, 55_000, (0, 1)),
    ("too_small", "DEL", 70_000, 80_000, (0, 1)),
    ("nothing", "DUP", 90_000, 95_000, (0, 1)),
    ("not_carrier", "DEL", 110_000, 115_000, (0, 1)),
]
CALLS = [
    ("DEL", 10_100, 15_000, "LowSQ", (0, 1), 1.2),   # weaker duplicate of the PASS call
    ("DEL", 10_200, 15_100, "PASS", (0, 1), 6.0),
    ("DUP", 30_050, 31_950, "LowSQ", (0, 1), 0.4),
    ("DUP", 50_000, 55_000, "PASS", (0, 1), 9.0),      # other type
    ("DEL", 70_000, 73_000, "PASS", (0, 1), 9.0),      # 30 % of the site
    ("DEL", 110_000, 115_000, "PASS", (0, 0), 9.0),    # non-carrier record
]


def _setup(tmp_path):
    _cohort(tmp_path / "cohort.vcf.gz", SITES)
    d = tmp_path / "samples" / "C"
    d.mkdir(parents=True)
    _calls(d / "C.filtered.vcf.gz", CALLS)
    (tmp_path / "ped.tsv").write_text("FAM\tF\t0\t0\t1\t1\nFAM\tM\t0\t0\t2\t1\nFAM\tC\tF\tM\t1\t2\n")
    return parse_ped(tmp_path / "ped.tsv")


def _rows(path):
    return {r["site_id"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def test_own_call_status(tmp_path):
    _setup(tmp_path)
    own = OwnCalls({"C": str(tmp_path / "samples" / "C" / "C.filtered.vcf.gz")})
    assert own.status("C", "chr1", 10_000, 15_000, "DEL") == "PASS"
    assert own.status("C", "chr1", 30_000, 32_000, "DUP") == "LowSQ"
    assert own.status("C", "chr1", 50_000, 55_000, "DEL") == "none"
    assert own.status("C", "chr1", 70_000, 80_000, "DEL") == "none"
    assert own.status("C", "chr1", 110_000, 115_000, "DEL") == "none"
    assert own.status("X", "chr1", 10_000, 15_000, "DEL") is None
    assert own.status("C", "chr2", 10_000, 15_000, "DEL") == "none"   # contig without calls
    own.close()


def test_genotype_only_label_and_counts(tmp_path):
    ped = _setup(tmp_path)
    [st] = run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", tmp_path / "s.tsv",
                            calls_dir=tmp_path / "samples")
    r = _rows(tmp_path / "m.tsv")
    assert {k: (v["inheritance"], v["child_call"]) for k, v in r.items()} == {
        "own_pass": ("denovo", "PASS"),
        "own_lowsq": ("denovo", "LowSQ"),
        "wrong_type": ("denovo_genotype_only", "none"),
        "too_small": ("denovo_genotype_only", "none"),
        "nothing": ("denovo_genotype_only", "none"),
        "not_carrier": ("denovo_genotype_only", "none"),
    }
    assert r["own_pass"]["dnq"] == "50" and r["nothing"]["dnq"] == ""
    assert (st.n_denovo, st.n_denovo_genotype_only) == (2, 4)
    summ = list(csv.DictReader(open(tmp_path / "s.tsv"), delimiter="\t"))
    assert (summ[0]["n_denovo"], summ[0]["n_denovo_genotype_only"]) == ("2", "4")


def test_low_confidence_stays_lowconf(tmp_path):
    """A de novo below the GQ floors keeps denovo_lowconf (with or without a call of
    the child): denovo + denovo_genotype_only is the old de novo count."""
    _cohort(tmp_path / "cohort.vcf.gz", SITES)
    ped = _setup(tmp_path)
    with pysam.VariantFile(str(tmp_path / "cohort.vcf.gz")) as vin:
        h = vin.header
        recs = list(vin)
    with pysam.VariantFile(str(tmp_path / "low.vcf.gz"), "wz", header=h) as out:
        for r in recs:
            r.samples["C"]["GQ"] = 25
            out.write(r)
    pysam.tabix_index(str(tmp_path / "low.vcf.gz"), preset="vcf", force=True)
    [st] = run_ped_annotate(tmp_path / "low.vcf.gz", ped, tmp_path / "m.tsv", None,
                            calls_dir=tmp_path / "samples")
    r = _rows(tmp_path / "m.tsv")
    assert {v["inheritance"] for v in r.values()} == {"denovo_lowconf"}
    assert r["own_pass"]["child_call"] == "PASS" and r["nothing"]["child_call"] == "none"
    assert (st.n_denovo, st.n_denovo_genotype_only) == (0, 0)


def test_without_calls_dir_unchanged(tmp_path):
    ped = _setup(tmp_path)
    [st] = run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", tmp_path / "s.tsv")
    r = _rows(tmp_path / "m.tsv")
    assert {v["inheritance"] for v in r.values()} == {"denovo"}
    assert {v["child_call"] for v in r.values()} == {""}
    assert (st.n_denovo, st.n_denovo_genotype_only) == (6, 0)
    summ = list(csv.DictReader(open(tmp_path / "s.tsv"), delimiter="\t"))
    assert summ[0]["n_denovo_genotype_only"] == "NA"
