import numpy as np

from svengers.breakpoints.insert_shift import dedupe_against, detect_insert_shift_dels
from svengers.config import BreakpointConfig
from svengers.extract.store import InsertModel
from svengers.variants import CNVCall

INSERT = InsertModel(mu=360.0, sigma=76.0, n_pairs=1_000_000)
BIN = 100


def _track(n_bins, rng, n_per_bin=15):
    n = np.full(n_bins, n_per_bin, dtype=np.uint32)
    means = rng.normal(INSERT.mu, INSERT.sigma / np.sqrt(n_per_bin), n_bins)
    return (means * n).astype(np.float64), n


def test_het_del_detected_with_size_estimate():
    rng = np.random.default_rng(0)
    s, n = _track(2000, rng)
    # het 160bp DEL centered at bin 500: spanning pairs shifted by +80 on average
    s[495:505] += 80 * n[495:505]
    calls = detect_insert_shift_dels("chr1", s, n, BIN, INSERT, BreakpointConfig())
    assert len(calls) == 1
    c = calls[0]
    assert c.svtype == "DEL" and c.source == "INSERT"
    assert 110 <= c.size <= 210          # ~2 x mean shift
    assert abs((c.start + c.end) // 2 - 50_000) <= 300
    assert c.gt == (0, 1) and c.gq > 0 and c.ins_pairs >= 50


def test_flat_track_no_calls():
    rng = np.random.default_rng(1)
    s, n = _track(2000, rng)
    assert detect_insert_shift_dels("chr1", s, n, BIN, INSERT, BreakpointConfig()) == []


def test_low_coverage_run_rejected():
    rng = np.random.default_rng(2)
    s, n = _track(2000, rng, n_per_bin=2)  # below min_n_per_bin
    s[500:510] += 200 * n[500:510]
    assert detect_insert_shift_dels("chr1", s, n, BIN, INSERT, BreakpointConfig()) == []


def test_dedupe_against_existing_del():
    ins = CNVCall("chr1", 49_950, 50_110, "DEL", "INSERT")
    junction_call = CNVCall("chr1", 49_960, 50_120, "DEL", "BOTH")
    unrelated = CNVCall("chr1", 200_000, 200_400, "DEL", "BOTH")
    assert dedupe_against([ins], [junction_call]) == []
    assert dedupe_against([ins], [unrelated]) == [ins]
