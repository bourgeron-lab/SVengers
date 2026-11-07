import numpy as np
import pytest

from vedette.config import GenotypeConfig
from vedette.genotype.likelihood import ContigData, Genotyper, expected_ratio
from vedette.variants import CNVCall


def _data(ratio, ploidy=2, baf_pos=None, baf_vaf=None, bin_size=1000, het_rate=None):
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
        het_rate=het_rate,
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


def test_paralog_dense_hets_do_not_refute_a_deletion():
    """24 hets in a 1.2 kb DEL (20 per kb, against 0.8 per kb on the contig)
    come from a paralog's reads: they must not veto a junction-supported DEL,
    while hets at the contig's own density still do."""
    rng = np.random.default_rng(5)
    ratio = rng.normal(1.0, 0.05, 100)
    ratio[50] = 0.7
    pos = np.arange(50_000, 51_200, 50)
    vaf = np.full(len(pos), 0.5, dtype=np.float32)
    gtr = Genotyper(GenotypeConfig())
    plain = _call(50_000, 51_200, "DEL", pe=8, sr=6)
    gtr.genotype(plain, _data(ratio, baf_pos=pos, baf_vaf=vaf))
    paralog = _call(50_000, 51_200, "DEL", pe=8, sr=6)
    gtr.genotype(paralog, _data(ratio, baf_pos=pos, baf_vaf=vaf, het_rate=8e-4))
    assert plain.gt == (0, 0)
    assert paralog.gt == (0, 1) and paralog.baf_n == 24
    off = Genotyper(GenotypeConfig(baf_paralog_p=0.0))
    call = _call(50_000, 51_200, "DEL", pe=8, sr=6)
    off.genotype(call, _data(ratio, baf_pos=pos, baf_vaf=vaf, het_rate=8e-4))
    assert call.gt == (0, 0)
    # 40 hets over 100 kb: the contig's density, still evidence against a DEL
    sparse = _call(0, 100_000, "DEL", pe=1, sr=0)
    wide = np.arange(10_000, 90_000, 2_000)
    gtr.genotype(sparse, _data(np.full(100, 0.7), baf_pos=wide,
                               baf_vaf=np.full(len(wide), 0.5, np.float32), het_rate=8e-4))
    assert sparse.gt == (0, 0)


def test_small_dup_below_bin_size_genotyped_from_junction():
    """A 100 bp DUP in 1 kb bins barely moves the depth; SR evidence must win."""
    rng = np.random.default_rng(6)
    ratio = rng.normal(1.0, 0.05, 100)  # flat depth (the DUP dilutes to ~5%)
    call = _call(50_000, 50_100, "DUP", pe=0, sr=8)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio))
    assert call.gt != (0, 0)


def test_overdispersion_flattens_het_hom_likelihood_gap():
    """The NB junction term must shrink the het-vs-hom log-likelihood gap a
    raw read count can produce (the measured Mendelian-violation mechanism),
    while barely touching the carrier-vs-noncarrier contrast."""
    from vedette.genotype.likelihood import _nb_loglik, _poisson_loglik

    lam_het, lam_hom, lam_bg = 2.1, 4.2, 0.175
    for count in (4, 8, 12):
        gap_pois = abs(_poisson_loglik(count, lam_hom) - _poisson_loglik(count, lam_het))
        gap_nb = abs(
            _nb_loglik(count, lam_hom, 0.75) - _nb_loglik(count, lam_het, 0.75)
        )
        assert gap_nb < gap_pois
        # carrier vs background stays strongly discriminative
        carrier_gap = _nb_loglik(count, lam_het, 0.75) - _nb_loglik(count, lam_bg, 0.75)
        assert carrier_gap > 3.0
    # alpha=0 recovers Poisson exactly
    assert _nb_loglik(5, 2.0, 0.0) == pytest.approx(_poisson_loglik(5, 2.0))


def test_haploid_del():
    rng = np.random.default_rng(5)
    ratio = rng.normal(0.5, 0.04, 100)
    ratio[20:40] = np.abs(rng.normal(0.01, 0.01, 20))
    call = _call(20_000, 40_000, "DEL", pe=10, sr=8)
    gtr = Genotyper(GenotypeConfig())
    gtr.genotype(call, _data(ratio, ploidy=1))
    assert call.gt == (1, 1)
    assert call.cn == 0
