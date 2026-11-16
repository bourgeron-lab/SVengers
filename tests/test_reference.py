import numpy as np

from svengers.reference import _contig_kind, _seq_tracks, bed_bin_overlap, build_bin_mask


def test_seq_tracks_gc_and_n():
    seq = b"GCGCGCGCGC" + b"ATATATATAT" + b"NNNNNNNNNN"
    gc, nfrac = _seq_tracks(seq, 10)
    assert np.allclose(gc, [1.0, 0.0, 0.0])
    assert np.allclose(nfrac, [0.0, 0.0, 1.0])


def test_seq_tracks_partial_last_bin_padded_as_n():
    gc, nfrac = _seq_tracks(b"GGGGG", 10)
    assert len(gc) == 1
    assert nfrac[0] == 0.5     # 5 padding bytes count as N
    assert gc[0] == 1.0        # GC fraction of the non-N half


def test_bed_bin_overlap():
    bed = {"chr1": np.array([[100, 250], [900, 1100]])}
    cover = bed_bin_overlap(bed, "chr1", n_bins=12, bin_size=100)
    assert cover[1] == 1.0                  # 100-200 fully covered
    assert np.isclose(cover[2], 0.5)        # 200-250
    assert cover[3] == 0.0
    assert cover[9] == 1.0 and np.isclose(cover[10], 1.0)
    assert bed_bin_overlap({}, "chr1", 5, 100).sum() == 0


def test_build_bin_mask():
    nfrac = np.array([0.0, 0.5, 0.0, 0.05], dtype=np.float32)
    bed = {"chr1": np.array([[200, 300]])}
    mask = build_bin_mask(4, 100, "chr1", nfrac, bed)
    assert mask.tolist() == [False, True, True, False]


def test_contig_kind():
    assert _contig_kind("chr7") == "autosome"
    assert _contig_kind("chrX") == "X"
    assert _contig_kind("chrY") == "Y"
    assert _contig_kind("chrM") == "MT"
    assert _contig_kind("chr1_KI270706v1_random") == "other"
    assert _contig_kind("HLA-A*01:01:01:01") == "other"


def test_interval_set_matches_across_contig_naming(tmp_path):
    """BED tracks mix UCSC and Ensembl naming; a 'chr'-less track must still
    answer for 'chr'-prefixed calls (the GHFC selfchain-LCR track is Ensembl-
    named, which silently zeroed every lcr_ovl until 2026-09-25)."""
    from svengers.filter.annotate import IntervalSet

    ensembl = tmp_path / "ensembl.bed"
    ensembl.write_text("1\t100\t200\n2\t0\t50\n")
    iv = IntervalSet.from_bed(ensembl)
    assert iv.overlap_fraction("chr1", 100, 200) == 1.0
    assert iv.overlap_fraction("1", 100, 200) == 1.0
    assert iv.overlap_fraction("chr1", 150, 250) == 0.5
    assert iv.overlap_fraction("chr3", 100, 200) == 0.0   # genuinely absent

    ucsc = tmp_path / "ucsc.bed"
    ucsc.write_text("chr1\t100\t200\n")
    iv2 = IntervalSet.from_bed(ucsc)
    assert iv2.overlap_fraction("1", 100, 200) == 1.0     # and the other way round
    assert iv2.overlap_fraction("chr1", 100, 200) == 1.0
