import pysam

from vedette.reference import Contig
from vedette.variants import CNVCall
from vedette.vcfout import clamp_end, write_vcf

CHRX_LEN = 156_040_895  # GRCh38 chrX; dataset_paper had a site ending 156,041,002


def test_clamp_end():
    assert clamp_end(1_000, 5_000, 10_000) == 5_000  # inside: untouched
    assert clamp_end(9_000, 11_000, 10_000) == 9_999  # written END = 10,000
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
    assert rec.start == 156_001_000
    assert rec.stop == CHRX_LEN  # END on the last base, not 156,041,001
    assert rec.info["SVLEN"] == -(CHRX_LEN - 1 - 156_001_000)
    assert "END=156040895;" in str(rec)
    # a call inside the contig keeps the usual END = end + 1 convention
    assert recs["DUP"].stop == 100_020_001
    assert recs["DUP"].info["SVLEN"] == 20_000
