import numpy as np

from vedette.config import GenotypeConfig
from vedette.genotype.likelihood import ContigData, Genotyper, expected_ratio
from vedette.variants import CNVCall


def _data(ratio, ploidy=2, baf_pos=None, baf_vaf=None, bin_size=1000):
    n = len(ratio)
    return ContigData(
        ratio=np.asarray(ratio, dtype=np.float64),
        mask=np.zeros(n, dtype=bool),
        baf_pos=np.asarray(baf_pos if baf_pos is not None else [], dtype=np.int64),
        baf_vaf=np.asarray(baf_vaf if baf_vaf is not None else [], dtype=np.float32),
        bin_size=bin_size,
        ploidy=ploidy,
        mean_depth=35.0,
        sigma0=0.1,
    )


def _call(start, end, svtype, pe=0, sr=0, source="BOTH"):
    return CNVCall(
        contig="chr1", start=start, end=end, svtype=svtype,
        source=source, pe_support=pe, sr_support=sr,
    )


def test_expected_ratio():
    assert expected_ratio(0, "DEL", 2) == 1.0
    assert expected_ratio(1, "DEL", 2) == 0.5
    assert expected_ratio(2, "DEL", 2) == 0.0
    assert expected_ratio(1, "DUP", 2) == 1.5
    assert expected_ratio(2, "DUP", 2) == 2.0
    assert expected_ratio(1, "DEL", 1) == 0.0


def test_het_del():
    rng = np.random.default_rng(0)
    ratio = rng.normal(1.0, 0.05, 100)
    ratio[40:60] = rng.normal(0.5, 0.05, 20)
    call = _call(40_000, 60_000, "DEL", pe=8, sr=6)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio))
    assert call.gt == (0, 1)
    assert call.gq >= 30
    assert call.cn == 1


def test_hom_ref_despite_junction_noise():
    """Normal depth + 1 stray supporting read must genotype 0/0."""
    rng = np.random.default_rng(2)
    ratio = rng.normal(1.0, 0.05, 100)
    call = _call(40_000, 60_000, "DUP", pe=1, sr=0)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio))
    assert call.gt == (0, 0)


def test_hom_del():
    rng = np.random.default_rng(3)
    ratio = rng.normal(1.0, 0.05, 100)
    ratio[40:60] = np.abs(rng.normal(0.01, 0.01, 20))
    call = _call(40_000, 60_000, "DEL", pe=20, sr=15)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio))
    assert call.gt == (1, 1)
    assert call.cn == 0


def test_baf_contradiction_lowers_het_del():
    """Clear 0.5 hets inside a putative het DEL contradict it."""
    rng = np.random.default_rng(4)
    ratio = rng.normal(0.70, 0.05, 100)  # ambiguous depth between CN1 and CN2
    pos = np.arange(10_000, 90_000, 2_000)
    vaf = np.full(len(pos), 0.5, dtype=np.float32)
    call_no_baf = _call(0, 100_000, "DEL", pe=1, sr=0)
    call_baf = _call(0, 100_000, "DEL", pe=1, sr=0)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call_no_baf, _data(ratio))
    gtr.genotype(call_baf, _data(ratio, baf_pos=pos, baf_vaf=vaf))
    # with 40 clear 0.5-VAF hets inside, the het-DEL hypothesis must lose ground
    assert call_baf.gt == (0, 0) or (
        call_no_baf.gt == (0, 1) and call_baf.gq < call_no_baf.gq
    )


def test_haploid_del():
    rng = np.random.default_rng(5)
    ratio = rng.normal(0.5, 0.04, 100)
    ratio[20:40] = np.abs(rng.normal(0.01, 0.01, 20))
    call = _call(20_000, 40_000, "DEL", pe=10, sr=8)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio, ploidy=1))
    assert call.gt == (1, 1)
    assert call.cn == 0
