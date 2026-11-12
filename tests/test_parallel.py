"""The cohort stages with jobs > 1 must write exactly what jobs = 1 writes."""

import csv
import gzip
import random

import numpy as np
import pysam
import pytest

CONTIGS = (("chr1", 5_000_000), ("chr2", 3_000_000), ("chr3", 2_000_000))
TRIOS = 3


def _samples():
    out = []
    for f in range(TRIOS):
        out += [f"F{f}P", f"F{f}M", f"F{f}C"]
    return out


def _pedigree(path):
    with open(path, "w") as fh:
        for f in range(TRIOS):
            fh.write(f"FAM{f}\tF{f}P\t0\t0\t1\t1\n")
            fh.write(f"FAM{f}\tF{f}M\t0\t0\t2\t1\n")
            fh.write(f"FAM{f}\tF{f}C\tF{f}P\tF{f}M\t1\t2\n")
    from vedette.ped import parse_ped

    return parse_ped(path)


def _filtered_vcf(path, sample, rng):
    header = pysam.VariantHeader()
    for name, length in CONTIGS:
        header.contigs.add(name, length=length)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("END", "1", "Integer", "End")
    header.filters.add("LowSQ", None, None, "below SQ threshold")
    header.filters.add("MaxN", None, None, "N overlap")
    for fid, num, typ in (("GT", "1", "String"), ("GQ", "1", "Integer"),
                          ("SQ", "1", "Float"), ("CN", "1", "Integer")):
        header.formats.add(fid, num, typ, fid)
    header.add_sample(sample)
    records = []
    for name, length in CONTIGS:
        # shared hotspots (recurrent sites, some straddling the 1 Mb regions)
        # plus private calls
        for centre in (400_000, 999_000, 1_000_500, 2_500_000):
            if centre < length and rng.random() < 0.7:
                start = centre + rng.randint(-300, 300)
                records.append((name, start, start + rng.choice((800, 5_000, 60_000))))
        for _ in range(12):
            start = rng.randrange(0, length - 200_000)
            records.append((name, start, start + rng.randint(200, 150_000)))
    order = {c: i for i, (c, _) in enumerate(CONTIGS)}
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, end in sorted(records, key=lambda r: (order[r[0]], r[1])):
            svtype = rng.choice(("DEL", "DEL", "DUP"))
            rec = out.new_record(contig=contig, start=start, stop=end,
                                 alleles=("N", f"<{svtype}>"))
            rec.info["SVTYPE"] = svtype
            flt = rng.choice(([], [], ["LowSQ"], ["LowSQ"], ["MaxN"]))
            for f in flt:
                rec.filter.add(f)
            fmt = rec.samples[sample]
            fmt["GT"] = rng.choice(((0, 1), (1, 1)))
            fmt["GQ"] = rng.randint(5, 99)
            fmt["SQ"] = round(rng.uniform(0.5, 9), 2)
            fmt["CN"] = rng.choice((0, 1)) if svtype == "DEL" else rng.choice((3, 4))
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _text(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.read()


@pytest.fixture
def cohort(tmp_path, monkeypatch):
    import vedette.parallel

    monkeypatch.setattr(vedette.parallel, "REGION_SIZE", 1_000_000)
    rng = random.Random(7)
    samples = _samples()
    paths = []
    for s in samples:
        d = tmp_path / "run" / "samples" / s
        d.mkdir(parents=True)
        p = d / f"{s}.filtered.vcf.gz"
        _filtered_vcf(p, s, rng)
        paths.append(p)
    return tmp_path, samples, paths, _pedigree(tmp_path / "ped.tsv")


def _genotype_tables(sites_vcf, samples, out_dir, rng):
    from vedette.cohort.genotyper import GENOTYPE_COLUMNS

    sites = [(r.id, r.contig, r.start, r.stop, r.info["SVTYPE"])
             for r in pysam.VariantFile(str(sites_vcf))]
    paths = []
    for s in samples:
        p = out_dir / f"{s}.gt.tsv"
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            for sid, contig, start, end, svtype in sites:
                gt = rng.choice(("0/0", "0/0", "0/1", "1/1", "./."))
                w.writerow(dict(zip(GENOTYPE_COLUMNS, [
                    sid, contig, start, end, svtype, s, gt, rng.randint(0, 99),
                    rng.choice(("", "-1", "1", "2", "3")), 0.9,
                    rng.choice(("", f"{rng.uniform(0.2, 1.8):.3f}")),
                    rng.randint(0, 20), rng.randint(0, 20), 0, 0, "",
                ], strict=False)))
        paths.append(p)
    return paths


def test_cohort_stages_parallel_equal_serial(cohort):
    from vedette.bench.trios import sweep
    from vedette.cohort.mendel import run_ped_annotate
    from vedette.cohort.merge import run_merge
    from vedette.cohort.metrics import run_cohort_filter
    from vedette.reference import Contig
    from vedette.report import _read_cohort_vcf, subthreshold_candidates

    tmp, samples, paths, ped = cohort
    contigs = [Contig(n, length, "autosome") for n, length in CONTIGS]

    n1 = run_merge(paths, contigs, tmp / "sites1.vcf.gz", jobs=1)
    n3 = run_merge(paths, contigs, tmp / "sites3.vcf.gz", jobs=3)
    assert n1 == n3 > 50
    assert _text(tmp / "sites1.vcf.gz") == _text(tmp / "sites3.vcf.gz")

    gts = _genotype_tables(tmp / "sites1.vcf.gz", samples, tmp, random.Random(11))
    c1 = run_cohort_filter(tmp / "sites1.vcf.gz", gts, tmp / "cohort1.vcf.gz",
                           pedigree=ped, jobs=1)
    c3 = run_cohort_filter(tmp / "sites1.vcf.gz", gts, tmp / "cohort3.vcf.gz",
                           pedigree=ped, jobs=3)
    assert c1 == c3
    assert _text(tmp / "cohort1.vcf.gz") == _text(tmp / "cohort3.vcf.gz")
    assert not list(tmp.glob(".cohort_parts_*"))   # temporary parts removed

    cohort_vcf = tmp / "cohort1.vcf.gz"
    s1 = run_ped_annotate(cohort_vcf, ped, tmp / "m1.tsv", tmp / "s1.tsv", jobs=1)
    s3 = run_ped_annotate(cohort_vcf, ped, tmp / "m3.tsv", tmp / "s3.tsv", jobs=3)
    assert [vars(x) for x in s1] == [vars(x) for x in s3]
    assert _text(tmp / "m1.tsv") == _text(tmp / "m3.tsv")
    assert _text(tmp / "s1.tsv") == _text(tmp / "s3.tsv")
    assert len(_text(tmp / "m1.tsv").splitlines()) > 50
    # the child's own calls: same rows serial and parallel, and some de novo of each kind
    own = tmp / "run" / "samples"
    o1 = run_ped_annotate(cohort_vcf, ped, tmp / "o1.tsv", tmp / "os1.tsv", jobs=1, calls_dir=own)
    o3 = run_ped_annotate(cohort_vcf, ped, tmp / "o3.tsv", tmp / "os3.tsv", jobs=3, calls_dir=own)
    assert [vars(x) for x in o1] == [vars(x) for x in o3]
    assert _text(tmp / "o1.tsv") == _text(tmp / "o3.tsv")
    inh = [ln.split("\t")[13] for ln in _text(tmp / "o1.tsv").splitlines()[1:]]
    assert "denovo_genotype_only" in inh
    assert sum(x.n_denovo + x.n_denovo_weak + x.n_denovo_genotype_only for x in o1) == sum(x.n_denovo for x in s1)

    assert sweep(cohort_vcf, ped, jobs=1) == sweep(cohort_vcf, ped, jobs=3)

    a = _read_cohort_vcf(cohort_vcf, with_cn=True, jobs=1)
    b = _read_cohort_vcf(cohort_vcf, with_cn=True, jobs=3)
    assert a[0] == b[0] and a[1] == b[1]
    for x, y in zip(a[2:], b[2:], strict=True):
        np.testing.assert_array_equal(x, y)

    run = tmp / "run"
    u, v = subthreshold_candidates(run, samples, jobs=1), subthreshold_candidates(run, samples, jobs=3)
    assert any(u.calls(s) for s in samples)
    assert all(u.calls(s) == v.calls(s) for s in samples)


def test_regions_follow_file_order_and_cover_each_record_once(cohort):
    from vedette.parallel import fetch_region, vcf_regions

    _tmp, _samples, paths, _ped = cohort
    path = paths[0]
    regions = vcf_regions(path)
    assert [r[0] for r in regions][:6] == ["chr1"] * 5 + ["chr2"]
    with pysam.VariantFile(str(path)) as vf:
        serial = [(r.contig, r.start, r.stop) for r in vf]
        split = [(r.contig, r.start, r.stop) for reg in regions for r in fetch_region(vf, reg)]
    assert split == serial
