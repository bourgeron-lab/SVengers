"""PASS calls the genotype prior made 0/0 report their best non-reference genotype.

The het prior (1e-4, ~9 nats) outweighs the few junction reads of a small
event on a 1 kb-depth genome: ~1,200 PASS calls per genome came out 0/0, most
of them true (HG002: 76 % of those < 300 bp, 94 % of 300 bp-1 kb), and then
dropped out of the cohort (merge loads carriers only) and of the family files.
"""

import gzip

import numpy as np
import pysam

from svengers.cohort.genotyper import StoreGenotyper, read_own_sites
from svengers.config import Config
from svengers.extract.store import EvidenceStore, EvidenceStoreWriter
from svengers.genotype.likelihood import ContigData, Genotyper
from svengers.reference import Contig
from svengers.variants import CNVCall
from svengers.vcfout import make_header, sv_record


def _flat(ploidy=2):
    n = 2000
    return ContigData(
        ratio=np.full(n, ploidy / 2.0), mask=np.zeros(n, bool),
        baf_pos=np.zeros(0, np.int64), baf_vaf=np.zeros(0), bin_size=1000,
        ploidy=ploidy, mean_depth=30.0, sigma0=0.1, het_rate=None,
    )


def _call(sr, pe=0, svtype="DEL", start=500_100, size=200):
    return CNVCall(contig="chr1", start=start, end=start + size, svtype=svtype,
                   source="JUNCTION", sr_support=sr, pe_support=pe)


def test_genotyper_reports_the_best_nonref_genotype():
    g = Genotyper(Config().genotype)
    weak = _call(sr=3)
    g.genotype(weak, _flat())
    assert weak.gt == (0, 0)                       # the prior wins
    assert weak.gt_alt == (0, 1) and weak.gq_alt > 0
    strong = _call(sr=6, pe=4)
    g.genotype(strong, _flat())
    assert strong.gt != (0, 0) and strong.gt_alt == strong.gt
    many = _call(sr=20, pe=15)                     # twice the het rate: homozygous
    g.genotype(many, _flat())
    assert many.gt_alt == (1, 1)
    male_x = _call(sr=2)
    g.genotype(male_x, _flat(ploidy=1))
    assert male_x.gt_alt == (1, 1)                 # haploid: the only non-reference genotype


def _store(path, sample="S1"):
    w = EvidenceStoreWriter(path, sample=sample, reference_id="x", bin_size=1000)
    w.set_meta(mean_depth=30.0, coverage=30.0, ploidy={"chr1": 2})
    depth = np.full(2000, 30.0, np.float32)
    w.write_depth("chr1", depth, depth, np.zeros(2000, bool))
    pos = np.arange(100_000, 1_000_000, 5_000, dtype=np.uint32)
    w.write_baf("chr1", pos, np.full(len(pos), 0.5, np.float32), np.full(len(pos), 30, np.uint16))
    return EvidenceStore(path)


def test_own_site_keeps_a_nonref_genotype(tmp_path):
    gtr = StoreGenotyper(_store(tmp_path / "S1.evidence.zarr"), Config())
    forced = gtr.genotype_site("s1", "chr1", 500_100, 500_300, "DEL")
    assert forced["gt"] == "0/0" and forced["own"] == 0
    own = gtr.genotype_site("s1", "chr1", 500_100, 500_300, "DEL", own=True)
    assert own["gt"] in ("0/1", "1/1") and own["own"] == 1
    gt, gq = gtr.nonref("chr1", 500_100, 500_300, "DEL", "JUNCTION", 0, 3, 0)
    assert gt == (0, 1) and 0 <= gq <= 99
    assert gtr.nonref("chr2", 500_100, 500_300, "DEL", "JUNCTION", 0, 3, 0) is None


def test_read_own_sites(tmp_path):
    p = tmp_path / "sites.members.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("svengers_site_1\tS1\nsvengers_site_1\tS10\nsvengers_site_2\tS10\n"
                 "svengers_site_3\tS1\n")
    assert read_own_sites(p, "S1") == {"svengers_site_1", "svengers_site_3"}
    assert read_own_sites(p, "S10") == {"svengers_site_1", "svengers_site_2"}


def _raw_calls(path, records):
    """records: (start, end, gt, sr), 0-based half-open"""
    hdr = make_header([Contig("chr1", 2_000_000, "autosome")], "S1")
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for i, (start, end, gt, sr) in enumerate(records):
            rec = sv_record(vf, "chr1", start, end, "DEL", f"svengers_DEL_{i}")
            rec.info["SVTYPE"], rec.info["SOURCE"] = "DEL", "JUNCTION"
            s = rec.samples["S1"]
            s["GT"], s["GQ"], s["PE"], s["SR"], s["CG"], s["CN"] = gt, 12, 0, sr, 0, 2
            vf.write(rec)


def test_filter_sets_pass_zero_zero_nonref(tmp_path):
    from svengers.filter.pipeline import run_filter

    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vcf, out = tmp_path / "S1.svengers.vcf.gz", tmp_path / "S1.filtered.vcf.gz"
    _raw_calls(vcf, [(500_100, 500_300, (0, 0), 3), (900_100, 900_400, (0, 1), 8)])
    counts = run_filter(vcf, store, out, Config())
    recs = list(pysam.VariantFile(str(out)))
    assert counts["pass"] == 2 and counts["nonref_gt"] == 1
    a, b = recs
    assert tuple(a.samples[0]["GT"]) == (0, 1) and a.info.get("NONREF_GT") is True
    assert tuple(b.samples[0]["GT"]) == (0, 1) and "NONREF_GT" not in b.info
    assert b.samples[0]["GQ"] == 12
    # switched off: the genotyper's call stands
    import dataclasses

    cfg = Config()
    off = dataclasses.replace(cfg, filter=dataclasses.replace(cfg.filter, pass_nonref_gt=False))
    out2 = tmp_path / "S1.off.vcf.gz"
    assert run_filter(vcf, store, out2, off)["nonref_gt"] == 0
    assert tuple(next(iter(pysam.VariantFile(str(out2)))).samples[0]["GT"]) == (0, 0)
