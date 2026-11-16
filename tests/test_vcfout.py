import pysam

from svengers.reference import Contig
from svengers.variants import CNVCall
from svengers.vcfout import clamp_end, span, spec_coords, write_vcf

CHRX_LEN = 156_040_895  # GRCh38 chrX; dataset_paper had a site ending 156,041,002


def test_clamp_end():
    assert clamp_end(1_000, 5_000, 10_000) == 5_000  # inside: untouched
    assert clamp_end(9_000, 11_000, 10_000) == 10_000  # written END = 10,000
    assert clamp_end(9_999, 11_000, 10_000) == 10_000  # never empty
    assert clamp_end(1_000, 5_000, None) == 5_000  # unknown contig: untouched


def test_depth_segment_last_bin_overhanging_contig(tmp_path):
    # 40 depth bins of 1 kb: the last bin [156,040,000, 156,041,000) runs
    # 105 bp past the end of chrX
    overhang = CNVCall(
        contig="chrX", start=156_001_000, end=156_041_000, svtype="DEL",
        source="DEPTH", n_bins=40, gt=(0, 1), gq=50, cn=1,
    )
    inside = CNVCall(
        contig="chrX", start=100_000_000, end=100_020_000, svtype="DUP",
        source="DEPTH", n_bins=20, gt=(0, 1), gq=50, cn=3,
    )
    out = tmp_path / "s.vcf.gz"
    write_vcf([overhang, inside], [Contig("chrX", CHRX_LEN, "X")], "S", out)

    recs = {r.info["SVTYPE"]: r for r in pysam.VariantFile(str(out))}
    rec = recs["DEL"]
    assert rec.pos == 156_001_000  # padding base before the first deleted one
    assert rec.stop == CHRX_LEN  # END on the last base, not 156,041,001
    assert rec.info["SVLEN"] == -(CHRX_LEN - 156_001_000)
    assert "END=156040895;" in str(rec)
    # inside the contig: POS = start, END = end (1-based), the VCF spec
    assert (recs["DUP"].pos, recs["DUP"].stop) == (100_000_000, 100_020_000)
    assert recs["DUP"].info["SVLEN"] == 20_000


def test_coordinates_round_trip(tmp_path):
    """POS is the padding base, END the last affected base, and span() gives
    back the internal half-open interval (genetics-viz issue: a junction
    deletion's POS sat one base past the breakpoint read from the clips)."""
    call = CNVCall(
        contig="chr1", start=39_537_353, end=39_543_232, svtype="DEL",
        source="JUNCTION", n_bins=6, gt=(0, 1), gq=50, cn=1,
    )
    out = tmp_path / "s.vcf.gz"
    write_vcf([call], [Contig("chr1", 248_956_422, "autosome")], "S", out)
    with pysam.VariantFile(str(out)) as vf:
        assert spec_coords(vf.header)
        (rec,) = vf
        assert (rec.pos, rec.stop, rec.info["SVLEN"]) == (39_537_353, 39_543_232, -5_879)
        assert span(rec, spec_coords(vf.header)) == (39_537_353, 39_543_232)


def test_legacy_vedette_coordinates(tmp_path):
    """A Vedette VCF written before v1.4.0 (POS = start + 1, END = end + 1)
    is read back to the same interval; foreign VCFs follow the spec."""
    header = pysam.VariantHeader()
    header.add_line("##source=vedette-1.3.3")
    header.contigs.add("chr1", length=1_000_000)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.info.add("SVTYPE", "1", "String", "x")
    header.info.add("SVLEN", "1", "Integer", "x")
    path = tmp_path / "old.vcf"
    with pysam.VariantFile(str(path), "w", header=header) as o:
        rec = o.new_record(contig="chr1", start=1_000, stop=1_100, alleles=("N", "<DEL>"))
        rec.info["SVTYPE"] = "DEL"
        rec.info["SVLEN"] = -100  # htslib: END = POS + 100 = 1,101
        o.write(rec)
    with pysam.VariantFile(str(path)) as vf:
        assert not spec_coords(vf.header)
        (rec,) = vf
        assert (rec.pos, rec.stop) == (1_001, 1_101)
        assert span(rec, False) == (1_000, 1_100)
    foreign = pysam.VariantHeader()
    foreign.add_line("##source=Manta")
    assert spec_coords(foreign)
