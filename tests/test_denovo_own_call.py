"""ped-annotate --calls-dir: a de novo without a call of the child at the site is genotype-only;
one with a non-PASS call of the child, a parental call or other carriers is weak."""

import csv

import pysam

from vedette.cohort.mendel import OwnCalls, run_ped_annotate
from vedette.ped import parse_ped

SAMPLES = ("C", "F", "M", "U")   # U: unrelated, carrier where the site says so


def _cohort(path, sites, unrelated=()):
    """sites: [(id, svtype, start, end, child_gt)], parents 0/0, GQ 50 everywhere; U 0/1 at
    the ``unrelated`` site ids."""
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
                r.samples[smp]["GT"] = (cgt if smp == "C" else
                                        (0, 1) if smp == "U" and sid in unrelated else (0, 0))
                r.samples[smp]["GQ"] = 50
            out.write(r)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _calls(path, calls, sample="C"):
    """calls: [(svtype, start, end, filter, gt, sq)] of one sample."""
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=1_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.filters.add("LowSQ", None, None, "low")
    h.formats.add("GT", "1", "String", "g")
    h.formats.add("SQ", "1", "Float", "s")
    h.add_sample(sample)
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        for t, s, e, flt, gt, sq in calls:
            r = out.new_record(contig="chr1", start=s - 1, stop=e, alleles=("N", f"<{t}>"))
            r.info["SVTYPE"] = t
            if flt != "PASS":
                r.filter.add(flt)
            r.samples[sample]["GT"] = gt
            r.samples[sample]["SQ"] = sq
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
    (tmp_path / "ped.tsv").write_text(
        "FAM\tF\t0\t0\t1\t1\nFAM\tM\t0\t0\t2\t1\nFAM\tC\tF\tM\t1\t2\nOTHER\tU\t0\t0\t1\t1\n"
    )
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
        "own_lowsq": ("denovo_weak", "LowSQ"),
        "wrong_type": ("denovo_genotype_only", "none"),
        "too_small": ("denovo_genotype_only", "none"),
        "nothing": ("denovo_genotype_only", "none"),
        "not_carrier": ("denovo_genotype_only", "none"),
    }
    assert r["own_pass"]["dnq"] == "50" and r["nothing"]["dnq"] == ""
    assert r["own_pass"]["parent_call"] == "NA,NA" and r["own_pass"]["other_carriers"] == "0"
    assert (st.n_denovo, st.n_denovo_weak, st.n_denovo_genotype_only) == (1, 1, 4)
    summ = list(csv.DictReader(open(tmp_path / "s.tsv"), delimiter="\t"))
    assert (summ[0]["n_denovo"], summ[0]["n_denovo_weak"], summ[0]["n_denovo_genotype_only"]) == (
        "1", "1", "4")


def test_parental_call_makes_it_weak(tmp_path):
    """A LowSQ call of the father at the site: the forced 0/0 missed him."""
    ped = _setup(tmp_path)
    d = tmp_path / "samples" / "F"
    d.mkdir()
    _calls(d / "F.filtered.vcf.gz", [("DEL", 10_050, 15_000, "LowSQ", (0, 1), 0.8)], sample="F")
    run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", None,
                     calls_dir=tmp_path / "samples")
    r = _rows(tmp_path / "m.tsv")["own_pass"]
    assert (r["inheritance"], r["parent_call"]) == ("denovo_weak", "LowSQ,NA")


def test_other_carriers_make_it_weak(tmp_path):
    """A carrier outside the family: allowed up to max_other_carriers."""
    ped = _setup(tmp_path)
    _cohort(tmp_path / "cohort.vcf.gz", SITES, unrelated={"own_pass"})
    for k, want in ((1, "denovo"), (0, "denovo_weak")):
        run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", None,
                         calls_dir=tmp_path / "samples", max_other_carriers=k)
        r = _rows(tmp_path / "m.tsv")["own_pass"]
        assert (r["inheritance"], r["other_carriers"]) == (want, "1")
    # the fraction of the pedigree's samples (4 here) wins when larger
    run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", None,
                     calls_dir=tmp_path / "samples", max_other_carriers=0, max_other_carrier_frac=0.25)
    assert _rows(tmp_path / "m.tsv")["own_pass"]["inheritance"] == "denovo"


def test_clingen_locus_lifts_the_carrier_limit(tmp_path):
    """A recurrent pathogenic de novo (PMS: terminal 22q13.3 over SHANK3 in half
    the children): carriers outside the family are allowed at a ClinGen locus of
    score 3 for the CNV's type; the call and parental-call rules still apply."""
    from vedette.clingen import ClinGen, Region

    ped = _setup(tmp_path)
    _cohort(tmp_path / "cohort.vcf.gz", SITES, unrelated={"own_pass", "own_lowsq"})

    def label(cg):
        run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", None,
                         calls_dir=tmp_path / "samples", max_other_carriers=0, clingen=cg)
        r = _rows(tmp_path / "m.tsv")
        return {k: (r[k]["inheritance"], r[k]["dosage_locus"]) for k in ("own_pass", "own_lowsq")}

    hi3 = ClinGen([], {}, [Region("SHANK3", "SHANK3", "chr1", 12_000, 13_000, 3, 1)])
    assert label(hi3) == {"own_pass": ("denovo", "SHANK3"), "own_lowsq": ("denovo_weak", "")}
    ts3 = ClinGen([], {}, [Region("G", "G", "chr1", 12_000, 13_000, 1, 3)])   # DEL needs HI 3
    assert label(ts3)["own_pass"] == ("denovo_weak", "")
    # the DUP site with a LowSQ call of the child stays weak at a TS 3 locus
    dup3 = ClinGen([], {}, [Region("G", "G", "chr1", 30_500, 31_000, 1, 3)])
    assert label(dup3)["own_lowsq"] == ("denovo_weak", "G")
    reg = ClinGen([Region("ISCA-9", "chr1 recurrent region", "chr1", 9_000, 30_000, 3, 0)], {})
    assert label(reg)["own_pass"] == ("denovo", "ISCA-9 chr1 recurrent region")
    assert label(None)["own_pass"] == ("denovo_weak", "")

    # a polymorphism in the gene: the carrier outside the family is a parent
    (tmp_path / "ped.tsv").write_text(
        "FAM\tF\t0\t0\t1\t1\nFAM\tM\t0\t0\t2\t1\nFAM\tC\tF\tM\t1\t2\n"
        "OTHER\tU\t0\t0\t1\t1\nOTHER\tV\t0\t0\t2\t1\nOTHER\tK\tU\tV\t1\t1\n")
    ped = parse_ped(tmp_path / "ped.tsv")
    assert label(hi3)["own_pass"] == ("denovo_weak", "")


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
    assert (st.n_denovo, st.n_denovo_weak, st.n_denovo_genotype_only) == (0, 0, 0)


def test_without_calls_dir_unchanged(tmp_path):
    ped = _setup(tmp_path)
    [st] = run_ped_annotate(tmp_path / "cohort.vcf.gz", ped, tmp_path / "m.tsv", tmp_path / "s.tsv")
    r = _rows(tmp_path / "m.tsv")
    assert {v["inheritance"] for v in r.values()} == {"denovo"}
    assert {v["child_call"] for v in r.values()} == {""}
    assert {v["parent_call"] for v in r.values()} == {""}
    assert (st.n_denovo, st.n_denovo_genotype_only) == (6, 0)
    summ = list(csv.DictReader(open(tmp_path / "s.tsv"), delimiter="\t"))
    assert summ[0]["n_denovo_genotype_only"] == summ[0]["n_denovo_weak"] == "NA"
