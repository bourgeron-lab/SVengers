"""Sex-chromosome karyotype from the evidence store and the PoN profiles."""

import numpy as np
import pytest
import zarr

from svengers.extract.store import EvidenceStore, EvidenceStoreWriter, _create_array
from svengers.karyotype import classify, is_anomaly, sex_chromosome_stats

BIN = 10_000
N = 3_000   # bins per contig (30 Mb)


def _pon(tmp_path):
    """chrX: bins 0-199 pseudo-autosomal (1.0 in both sexes), the rest 1.0 in
    women / 0.5 in men; chrY: bins 0-999 a repeat seen in women too."""
    root = zarr.open(str(tmp_path / "pon.zarr"), mode="w")
    for c in ("chrX", "chrY"):
        g = root.create_group(c)
        if c == "chrX":
            f = np.ones(N)
            m = np.full(N, 0.5)
            m[:200] = 1.0
        else:
            f = np.zeros(N)
            f[:1000] = 0.4
            m = np.full(N, 0.5)
        for name, arr in (("median_female", f), ("median_male", m),
                          ("mad_female", np.full(N, 0.05)), ("mad_male", np.full(N, 0.05))):
            _create_array(g, name, arr, N)
    return root


def _store(tmp_path, x_copies, y_copies, het_x=None, x_baf_shift=0.0, seed=0, auto=None):
    rng = np.random.default_rng(seed)
    path = tmp_path / f"S{x_copies}{y_copies}.evidence.zarr"
    w = EvidenceStoreWriter(path, sample="S", reference_id="x", bin_size=BIN)
    # auto: {autosome: copies} beyond chr2/7/12 (a GC bias, a trisomy)
    contigs = [("chr2", 2), ("chr7", 2), ("chr12", 2), *(auto or {}).items()]
    for c, cn in (*contigs, ("chrX", x_copies), ("chrY", y_copies)):
        depth = (15.0 * cn + rng.normal(0, 1.0, N)).clip(0).astype(np.float32)
        if c == "chrX":
            depth[:200] = 30.0                      # pseudo-autosomal: 2 copies in both sexes
        w.write_depth(c, depth, depth, np.zeros(N, dtype=bool))
        n = 20_000
        pos = np.sort(rng.integers(0, N * BIN, n)).astype(np.uint32)
        het_rate = 0.6 if c != "chrX" else (het_x if het_x is not None else (0.6 if x_copies >= 2 else 0.003))
        het = rng.random(n) < het_rate
        side = rng.choice([-1, 1], n)
        shift = x_baf_shift if c == "chrX" else 0.0
        vaf = np.where(het, 0.5 + side * shift + rng.normal(0, 0.07, n), 1.0)
        w.write_baf(c, pos, np.clip(vaf, 0.0, 1.0).astype(np.float32), np.zeros(n, np.uint16))
    return EvidenceStore(path)


@pytest.mark.parametrize("x,y,expected", [(2, 0, "XX"), (1, 1, "XY"), (1, 0, "X0"),
                                          (2, 1, "XXY"), (1, 2, "XYY")])
def test_karyotypes(tmp_path, x, y, expected):
    st = sex_chromosome_stats(_store(tmp_path, x, y), _pon(tmp_path))
    assert st["x_cn"] == pytest.approx(x, abs=0.1)
    assert (st["y_cn"] or 0) == pytest.approx(y, abs=0.1)
    k = classify(st)
    assert k["karyotype"] == expected
    assert is_anomaly(k) == (expected not in ("XX", "XY"))


def test_triple_x_and_its_baf(tmp_path):
    st = sex_chromosome_stats(_store(tmp_path, 3, 0, x_baf_shift=1 / 6), _pon(tmp_path))
    assert classify(st)["karyotype"] == "XXX"
    assert st["x_baf"] == pytest.approx(1 / 6, abs=0.03)


def test_gc_bias_voids_an_abnormal_reading(tmp_path):
    pon = _pon(tmp_path)
    # GC-rich autosomes at 3 copies, AT-rich at 1.6, chrX read at 1.6 in a woman
    biased = {"chr4": 1.6, "chr5": 1.6, "chr13": 1.6, "chr18": 1.6,
              "chr16": 3, "chr17": 3, "chr19": 3, "chr22": 3}
    st = sex_chromosome_stats(_store(tmp_path, 1.6, 0, het_x=0.6, auto=biased), pon)
    assert st["gc_contrast"] == pytest.approx(1.4, abs=0.1)
    assert classify(st, "F")["karyotype"] == "?" and "biais GC" in classify(st, "F")["note"]
    # an XY reading survives the same bias
    st = sex_chromosome_stats(_store(tmp_path, 1, 1, auto=biased), pon)
    assert classify(st, "M") == {"karyotype": "XY", "note": ""}
    # a trisomy 18 is not a GC bias
    tri = {c: 2 for c in biased} | {"chr18": 3}
    st = sex_chromosome_stats(_store(tmp_path, 2, 1, seed=1, auto=tri), pon)
    assert st["gc_contrast"] == pytest.approx(0, abs=0.1)
    assert classify(st, "M")["karyotype"] == "XXY"


def test_pedigree_discordance_contamination_and_mosaic():
    xy = {"x_cn": 1.0, "y_cn": 0.98, "x_het": 0.01, "x_baf": 0.38}
    assert "pedigree" in classify(xy, "F")["note"]
    assert classify(xy, "M") == {"karyotype": "XY", "note": ""}
    assert "contamination" in classify({**xy, "x_het": 0.4}, "M")["note"]
    assert classify({"x_cn": 1.7, "y_cn": 0.0, "x_het": 1.0, "x_baf": 0.06})["karyotype"] == "intermédiaire"
    assert "perte partielle du Y" in classify({**xy, "y_cn": 0.55})["note"]


@pytest.mark.parametrize("kar,expected", [
    ("XX", set()), ("XY", set()), ("?", set()), ("intermédiaire", set()),
    ("XXY", {"X"}), ("XXX", {"X"}), ("X0", {"X"}), ("XYY", {"Y"}), ("XXYY", {"X", "Y"}),
])
def test_aneuploid_sex_chromosomes(kar, expected):
    from svengers.karyotype import aneuploid_sex_chromosomes

    assert aneuploid_sex_chromosomes(kar) == expected


@pytest.mark.parametrize("kar,expected", [
    ("XX", {"chrX": 2, "chrY": 0}), ("XY", {"chrX": 1, "chrY": 1}), ("X0", {"chrX": 1, "chrY": 0}),
    ("XXY", {"chrX": 2, "chrY": 1}), ("XXX", {"chrX": 3, "chrY": 0}), ("XYY", {"chrX": 1, "chrY": 2}),
    ("?", None), ("intermédiaire", None), ("", None),
])
def test_sex_ploidy(kar, expected):
    from svengers.karyotype import sex_ploidy

    assert sex_ploidy(kar) == expected


@pytest.mark.parametrize("kind,ploidy,expected", [
    ("autosome", 2, ("all", 1.0)), ("X", 1, ("male", 0.5)), ("X", 2, ("female", 1.0)),
    ("X", 3, ("female", 1.0)), ("Y", 1, ("male", 0.5)), ("Y", 2, ("male", 0.5)), ("Y", 0, ("female", 0.0)),
])
def test_pon_profile_level_is_the_profile_ploidy(kind, ploidy, expected):
    from svengers.call import pon_profile

    assert pon_profile(kind, ploidy) == expected


def test_pon_normalize_keeps_an_aneuploid_on_its_own_baseline():
    """An XXX chrX (ratio 1.5) on the female profile reads 1.5 = ploidy 3 / 2;
    an XYY chrY (ratio 1.0) on the male profile reads 1.0 = ploidy 2 / 2."""
    from svengers.call import pon_profile
    from svengers.cohort.pon import PanelOfNormals

    class P(PanelOfNormals):
        def __init__(self, med):
            self.med = med

        def profile(self, contig, sex="all"):
            return np.full(100, self.med[sex]), np.full(100, 0.02)

    pon = P({"female": 1.0, "male": 0.5})
    mask = np.zeros(100, bool)
    for kind, ploidy, ratio in (("X", 3, 1.5), ("X", 2, 1.0), ("X", 1, 0.5), ("Y", 2, 1.0), ("Y", 1, 0.5)):
        sex, level = pon_profile(kind, ploidy)
        out, _ = pon.normalize("c", np.full(100, ratio), mask, sex=sex, level=level, denoise=False)
        assert out.mean() == pytest.approx(ploidy / 2)


def test_karyotype_from_call_arrays_matches_the_store(tmp_path):
    from svengers.karyotype import ArraySource

    pon = _pon(tmp_path)
    st = _store(tmp_path, 2, 1)
    contigs = ["chr2", "chr7", "chr12", "chrX", "chrY"]
    raw = {c: np.asarray(st.depth(c)[0]) for c in contigs}
    mask = {c: np.asarray(st.depth(c)[2]) for c in contigs}
    baf = {}
    for c in contigs:
        b = st.baf(c, 0, 1 << 30)
        baf[c] = (np.asarray(b["pos"]), np.asarray(b["vaf"]), np.zeros(len(b["pos"]), np.uint16))
    a = sex_chromosome_stats(ArraySource("S", BIN, raw, mask, baf), pon)
    b = sex_chromosome_stats(st, pon)
    assert a == b
    assert classify(a)["karyotype"] == "XXY"


def test_pon_keeps_sex_aneuploidies_out_of_the_sex_profiles():
    from svengers.cohort.pon import _sex_of_store

    class S:
        def __init__(self, x, y=None):
            self.meta = {"ploidy": {"chrX": x, **({} if y is None else {"chrY": y})}}

    assert _sex_of_store(S(1, 1)) == "male" and _sex_of_store(S(1)) == "male"
    assert _sex_of_store(S(2, 0)) == "female" and _sex_of_store(S(2)) == "female"
    assert _sex_of_store(S(2, 1)) is None      # XXY
    assert _sex_of_store(S(3, 0)) is None      # XXX
    assert _sex_of_store(S(1, 2)) is None      # XYY
    assert _sex_of_store(S(1, 0)) is None      # X0


def test_vcf_header_carries_the_karyotype(tmp_path):
    import pysam

    from svengers.reference import Contig
    from svengers.vcfout import write_vcf

    out = tmp_path / "s.vcf.gz"
    write_vcf([], [Contig("chrX", 1_000_000, "X")], "S", out, karyotype="XXY")
    assert "##svengers_karyotype=XXY" in str(pysam.VariantFile(str(out)).header)


def _autosomes(levels, baf_dev=None, het=None, n_bins=40_000, seed=0):
    """ratio/mask/baf dicts at 1 kb bins; levels: {contig: copies}."""
    rng = np.random.default_rng(seed)
    ratio, mask, baf = {}, {}, {}
    for c, k in levels.items():
        ratio[c] = (k / 2 + rng.normal(0, 0.08, n_bins)).astype(np.float32)
        mask[c] = np.zeros(n_bins, bool)
        n = 30_000
        pos = np.sort(rng.integers(0, n_bins * 1000, n)).astype(np.uint32)
        rate = (het or {}).get(c, 0.6)
        is_het = rng.random(n) < rate
        dev = (baf_dev or {}).get(c, 0.0)
        vaf = np.where(is_het, 0.5 + rng.choice([-1, 1], n) * dev + rng.normal(0, 0.06, n), 1.0)
        baf[c] = (pos, np.clip(vaf, 0, 1).astype(np.float32), np.zeros(n, np.uint16))
    return ratio, mask, baf


def test_whole_chromosome_trisomy_and_monosomy():
    from svengers.karyotype import autosomal_aneuploidies

    lv = {"chr1": 2, "chr2": 2, "chr12": 3, "chr21": 3, "chr22": 1, "chr7": 2}
    r, m, b = _autosomes(lv, baf_dev={"chr12": 1 / 6, "chr21": 1 / 6}, het={"chr22": 0.003})
    assert autosomal_aneuploidies(r, m, b, 1000, list(lv)) == {"chr12": 3, "chr21": 3, "chr22": 1}


def test_trisomy_of_a_small_acrocentric_with_unmasked_gaps():
    """chr21 as the call sees it before the PoN: ~37 informative windows, the
    p-arm flank unmasked but uncovered (CN 0) or erratic (0.6, 4.9)."""
    from svengers.karyotype import autosomal_aneuploidies

    lv = {"chr1": 2, "chr2": 2, "chr21": 3}
    r, m, b = _autosomes(lv, baf_dev={"chr21": 1 / 6}, n_bins=37_000)
    r["chr21"][:2_000] = 0.0                       # two uncovered windows
    r["chr21"][2_000:3_000] = 0.31                 # erratic flank
    r["chr21"][3_000:4_000] = 2.43
    assert autosomal_aneuploidies(r, m, b, 1000, list(lv)) == {"chr21": 3}


def test_depth_alone_is_not_an_aneuploidy():
    """A chromosome read at 3 copies without the BAF of three (a wave, a GC
    bias) or a gain over part of a chromosome keeps the diploid baseline."""
    from svengers.karyotype import autosomal_aneuploidies

    lv = {"chr1": 2, "chr19": 3, "chr18": 2}
    r, m, b = _autosomes(lv)                       # chr19 BAF stays at 0.5
    r["chr18"][: 15_000] = 1.5                     # 15 of 40 Mb gained
    assert autosomal_aneuploidies(r, m, b, 1000, list(lv)) == {}
    # a mosaic level (2.5 copies) is not an integer baseline either
    r, m, b = _autosomes({"chr1": 2, "chr12": 2.5}, baf_dev={"chr12": 0.08})
    assert autosomal_aneuploidies(r, m, b, 1000, ["chr1", "chr12"]) == {}


def test_vcf_header_carries_the_aneuploidies(tmp_path):
    import pysam

    from svengers.reference import Contig
    from svengers.vcfout import write_vcf

    out = tmp_path / "s.vcf.gz"
    write_vcf([], [Contig("chr21", 1_000_000, "autosome")], "S", out, aneuploidies={"chr21": 3})
    assert "##svengers_aneuploidy=chr21:3" in str(pysam.VariantFile(str(out)).header)
