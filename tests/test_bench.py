"""Integration tests of the bench harness on tiny synthetic data (real truvari)."""

import shutil
import sys
from pathlib import Path

import pysam
import pytest

from vedette.bench.compare import run_compare, wisecondorx_to_vcf
from vedette.bench.plots import plot_size_stratified, plot_trio_sweep
from vedette.bench.truth import Stratifier, _size_class, evaluate
from vedette.reference import Contig

truvari_missing = (
    shutil.which("truvari") is None
    and not Path(sys.executable).with_name("truvari").exists()
)

CONTIGS = [Contig("chr1", 1_000_000, "autosome")]


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


@pytest.fixture()
def truth_and_bed(tmp_path):
    truth = tmp_path / "truth.vcf.gz"
    _write(truth, _header(), [
        # 500 bp DEL (300bp-1kb class), 5 kb DEL (1-10kb class)
        ("chr1", 100_000, 100_501, ("A" + "T" * 500, "A"), "DEL", -500, "t1"),
        ("chr1", 300_000, 305_001, ("A" + "T" * 5000, "A"), "DEL", -5000, "t2"),
    ])
    bed = tmp_path / "conf.bed"
    bed.write_text("chr1\t0\t1000000\n")
    return truth, bed


def test_size_class_boundaries():
    assert _size_class(50) == "50-300bp"
    assert _size_class(999) == "300bp-1kb"
    assert _size_class(10_000) == "10-100kb"
    assert _size_class(2_000_000) == ">1Mb"


@pytest.mark.skipif(truvari_missing, reason="truvari not installed")
def test_evaluate_stratified(tmp_path, truth_and_bed):
    truth, bed = truth_and_bed
    calls = tmp_path / "calls.vcf.gz"
    _write(calls, _header(with_symbolic=True), [
        ("chr1", 100_005, 100_500, ("N", "<DEL>"), "DEL", -495, "c1"),   # TP (500bp)
        ("chr1", 600_000, 600_400, ("N", "<DEL>"), "DEL", -400, "c_fp"),  # FP
    ])
    rows = evaluate(calls, truth, bed, tmp_path / "eval", stratifier=Stratifier())
    by_key = {(r["svtype"], r["stratum_kind"], r["stratum"]): r for r in rows}

    overall = by_key[("DEL", "all", "all")]
    assert overall["tp_comp"] == 1 and overall["fp"] == 1
    assert overall["tp_base"] == 1 and overall["fn"] == 1  # 5kb DEL missed
    assert float(overall["recall"]) == 0.5 and float(overall["precision"]) == 0.5

    assert by_key[("DEL", "size", "1-10kb")]["fn"] == 1        # the missed 5kb DEL
    assert by_key[("DEL", "size", "300bp-1kb")]["tp_base"] == 1
    assert by_key[("DEL", "context", "clean")]["tp_comp"] == 1

    plot_size_stratified(rows, tmp_path / "f1.png")
    assert (tmp_path / "f1.png").stat().st_size > 0


@pytest.mark.skipif(truvari_missing, reason="truvari not installed")
def test_compare_with_wisecondorx_entry(tmp_path, truth_and_bed):
    truth, bed = truth_and_bed
    wcx = tmp_path / "aberrations.bed"
    wcx.write_text(
        "chr\tstart\tend\tratio\tzscore\ttype\tsource\n"
        "1\t99900\t100600\t-1.0\t-8.0\tprimary\tloss\n"     # matches t1
        "1\t700000\t710000\t0.58\t9.0\tprimary\tgain\n"     # FP DUP
    )
    calls = tmp_path / "calls.vcf.gz"
    _write(calls, _header(with_symbolic=True), [
        ("chr1", 100_005, 100_500, ("N", "<DEL>"), "DEL", -495, "c1"),
        ("chr1", 299_990, 305_000, ("N", "<DEL>"), "DEL", -5010, "c2"),
    ])
    rows = run_compare(
        [("vedette", "vcf", str(calls)), ("wisecondorx", "wisecondorx", str(wcx))],
        truth, bed, CONTIGS, tmp_path / "cmp",
    )
    by = {
        (r["caller"], r["svtype"], r["stratum"]): r
        for r in rows if r["stratum_kind"] == "all"
    }
    assert float(by[("vedette", "DEL", "all")]["recall"]) == 1.0
    assert float(by[("wisecondorx", "DEL", "all")]["recall"]) == 0.5
    assert (tmp_path / "cmp" / "comparison.tsv").exists()


def test_wisecondorx_conversion_adds_chr_prefix(tmp_path):
    bed = tmp_path / "ab.bed"
    bed.write_text(
        "chr\tstart\tend\tratio\tzscore\ttype\tsource\n"
        "5\t1000\t2000\t-1.0\t-5.0\tprimary\tloss\n"
    )
    out = tmp_path / "out.vcf.gz"
    n = wisecondorx_to_vcf(bed, out, [Contig("chr5", 100_000, "autosome")])
    assert n == 1
    with pysam.VariantFile(str(out)) as vcf:
        rec = next(iter(vcf))
        assert rec.contig == "chr5" and rec.info["SVTYPE"] == "DEL"


def test_plot_trio_sweep(tmp_path):
    rows = [
        {"threshold": t, "field": "GQ", "mvr": 0.01 / (t + 1), "transmission": 0.6 - 0.01 * t}
        for t in (0, 10, 20)
    ]
    plot_trio_sweep(rows, tmp_path / "sweep.png")
    assert (tmp_path / "sweep.png").stat().st_size > 0