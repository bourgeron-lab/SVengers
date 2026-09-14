import numpy as np

from vedette.qc import WAVE_MAX_K, contig_wave_score, sample_wave_score

BIN = 1000
N = 60_000  # 60 Mb contig at 1 kb bins


def _flat(rng, n=N):
    return rng.normal(1.0, 0.05, n)


def _wavy(rng, n=N, period_bins=800, amp=0.15):
    x = np.arange(n)
    return 1.0 + amp * np.sin(2 * np.pi * x / period_bins) + rng.normal(0, 0.05, n)


def test_flat_profile_scores_low():
    rng = np.random.default_rng(0)
    s = contig_wave_score(_flat(rng), np.zeros(N, bool), BIN)
    # white noise: energy fraction ~= band width / spectrum width ~ 1%
    assert s is not None and s < 0.10


def test_wavy_profile_scores_high():
    rng = np.random.default_rng(0)
    s = contig_wave_score(_wavy(rng), np.zeros(N, bool), BIN)  # 800 kb period
    assert s is not None and s > 0.5


def test_wave_score_orders_amplitudes():
    rng = np.random.default_rng(1)
    weak = contig_wave_score(_wavy(rng, amp=0.05), np.zeros(N, bool), BIN)
    strong = contig_wave_score(_wavy(rng, amp=0.3), np.zeros(N, bool), BIN)
    assert weak < strong


def test_masked_stretches_do_not_fabricate_waves():
    rng = np.random.default_rng(2)
    ratio = _flat(rng)
    mask = np.zeros(N, bool)
    mask[::7] = True          # scattered masked bins
    mask[10_000:12_000] = True  # one long masked block (interp would ramp here)
    s = contig_wave_score(ratio, mask, BIN)
    assert s is not None and s < 0.10


def test_short_contig_returns_none():
    rng = np.random.default_rng(3)
    assert contig_wave_score(_flat(rng, 3_000), np.zeros(3_000, bool), BIN) is None


def test_cn_step_is_not_a_wave():
    # a real 5 Mb heterozygous deletion is a step, not periodicity: the linear
    # detrend won't remove it entirely but it must not reach the wavy range
    rng = np.random.default_rng(4)
    ratio = _flat(rng)
    ratio[20_000:25_000] *= 0.5
    s = contig_wave_score(ratio, np.zeros(N, bool), BIN)
    wavy = contig_wave_score(_wavy(rng), np.zeros(N, bool), BIN)
    assert s < wavy


def test_sample_score_aggregates_and_gates_on_contig_count():
    rng = np.random.default_rng(5)
    ratio = {f"chr{i}": _wavy(rng) for i in range(1, 13)}
    mask = {k: np.zeros(N, bool) for k in ratio}
    score, per_contig = sample_wave_score(ratio, mask, BIN)
    assert score is not None and score > 0.5 and len(per_contig) == 12
    one = {"chr20": ratio["chr1"]}
    score_one, per_one = sample_wave_score(one, {"chr20": mask["chr1"]}, BIN)
    assert score_one is None and "chr20" in per_one


def test_band_definition_sane():
    assert WAVE_MAX_K == 50
