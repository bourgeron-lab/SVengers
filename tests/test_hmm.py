import numpy as np
import pytest

from vedette.segment import hmm


def _ratio(n, rng, sigma=0.08):
    return rng.normal(1.0, sigma, n).astype(np.float64)


def test_het_del_and_dup_detected():
    rng = np.random.default_rng(42)
    ratio = _ratio(2000, rng)
    ratio[500:600] = rng.normal(0.5, 0.08, 100)   # het DEL
    ratio[1500:1550] = rng.normal(1.5, 0.08, 50)  # het DUP
    mask = np.zeros(2000, dtype=bool)
    segments, path = hmm.segment_contig(ratio, mask, ploidy=2, max_cn=5, transition=1e-6)
    kinds = {(s.cn, s.start_bin, s.end_bin) for s in segments}
    assert any(cn == 1 and abs(s - 500) <= 2 and abs(e - 600) <= 2 for cn, s, e in kinds)
    assert any(cn == 3 and abs(s - 1500) <= 2 and abs(e - 1550) <= 2 for cn, s, e in kinds)
    assert all(s.posterior > 0.9 for s in segments)


def test_hom_del_zero_depth():
    rng = np.random.default_rng(1)
    ratio = _ratio(1000, rng)
    ratio[300:350] = np.abs(rng.normal(0.01, 0.01, 50))
    mask = np.zeros(1000, dtype=bool)
    segments, _ = hmm.segment_contig(ratio, mask, ploidy=2, max_cn=5, transition=1e-6)
    assert any(s.cn == 0 and abs(s.start_bin - 300) <= 2 for s in segments)


def test_masked_gap_does_not_split_segment():
    rng = np.random.default_rng(7)
    ratio = _ratio(1000, rng)
    ratio[400:500] = rng.normal(0.5, 0.08, 100)
    mask = np.zeros(1000, dtype=bool)
    mask[440:460] = True          # masked gap inside the DEL
    ratio[440:460] = 5.0          # garbage values under the mask must be ignored
    segments, _ = hmm.segment_contig(ratio, mask, ploidy=2, max_cn=5, transition=1e-6)
    dels = [s for s in segments if s.cn == 1]
    assert len(dels) == 1
    assert dels[0].start_bin <= 402 and dels[0].end_bin >= 498


def test_flat_genome_no_calls():
    rng = np.random.default_rng(3)
    ratio = _ratio(5000, rng)
    mask = np.zeros(5000, dtype=bool)
    segments, _ = hmm.segment_contig(ratio, mask, ploidy=2, max_cn=5, transition=1e-8)
    assert segments == []


def test_haploid_baseline():
    rng = np.random.default_rng(11)
    ratio = rng.normal(0.5, 0.06, 1000)  # male chrX baseline
    ratio[200:260] = np.abs(rng.normal(0.01, 0.01, 60))  # hemizygous DEL
    mask = np.zeros(1000, dtype=bool)
    segments, _ = hmm.segment_contig(ratio, mask, ploidy=1, max_cn=5, transition=1e-6)
    assert any(s.cn == 0 and abs(s.start_bin - 200) <= 2 for s in segments)
    # baseline CN1 must not be reported as a genome-wide DEL
    assert not any(s.cn == 1 for s in segments)


def test_baf_emissions_shift_ambiguous_dup():
    """Clear 1/3-2/3 het VAFs must push an ambiguous depth ratio toward CN3."""
    rng = np.random.default_rng(5)
    n = 400
    ratio = _ratio(n, rng)
    ratio[100:200] = rng.normal(1.30, 0.10, 100)  # between CN2 (1.0) and CN3 (1.5)
    mask = np.zeros(n, dtype=bool)
    pos = np.arange(100_000, 200_000, 500, dtype=np.uint32)
    vaf = np.where(np.arange(len(pos)) % 2 == 0, 0.33, 0.67).astype(np.float32)
    baf_emit = hmm.baf_log_emissions(pos, vaf, n, 1000, max_cn=5)
    without, _ = hmm.segment_contig(ratio, mask, 2, 5, 1e-6)
    with_baf, _ = hmm.segment_contig(ratio, mask, 2, 5, 1e-6, baf_emit=baf_emit)
    assert any(s.cn == 3 and s.start_bin <= 102 for s in with_baf)


def test_fit_sigma0_hemizygous_contig():
    """A male chrX (baseline 0.5) gives the same CN2-unit sigma0 as an autosome
    of the same sample, and a single noisy bin is not a hemizygous DEL."""
    rng = np.random.default_rng(4)
    sd2 = 0.08
    auto = rng.normal(1.0, sd2, 50_000)
    chrx = rng.normal(0.5, sd2 * np.sqrt(0.5), 50_000)
    mask = np.zeros(50_000, dtype=bool)
    s_auto = hmm.fit_sigma0(auto, mask)
    s_x = hmm.fit_sigma0(chrx, mask, ploidy=1)
    assert s_auto == pytest.approx(sd2, rel=0.05)
    assert s_x == pytest.approx(sd2, rel=0.05)
    assert hmm.fit_sigma0(chrx, mask, ploidy=0) == 0.15
    # no CN change on a hemizygous contig of pure noise
    segments, path = hmm.segment_contig(chrx, mask, ploidy=1, max_cn=5, transition=1e-4)
    assert not segments and set(np.unique(path)) == {1}
