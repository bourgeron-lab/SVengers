import numpy as np
import pytest

from vedette.cohort.merge import MemberCall, SiteCluster, _compatible, cluster_calls
from vedette.cohort.metrics import hwe_exact


def _mc(sample, start, end, precise=True, gq=99):
    return MemberCall(sample=sample, start=start, end=end, precise=precise, gq=gq)


def test_identical_calls_cluster():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 1010, 5005), "DEL")]
    clusters = cluster_calls(calls)
    assert len(clusters) == 1
    assert clusters[0].members[0].sample == "A"
    assert len(clusters[0].members) == 2


def test_del_dup_never_merge():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 1000, 5000), "DUP")]
    assert len(cluster_calls(calls)) == 2


def test_distant_calls_stay_separate():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 50_000, 55_000), "DEL")]
    assert len(cluster_calls(calls)) == 2


def test_no_transitive_chaining():
    """A drifting chain must not fuse into one mega-cluster (AF inflation guard)."""
    calls = []
    start, end = 10_000, 20_000
    for i in range(30):
        calls.append((_mc(f"S{i}", start, end, precise=False), "DEL"))
        start += 1_500  # each next call shifts; pairwise-adjacent are compatible
        end += 1_500
    clusters = cluster_calls(sorted(calls, key=lambda t: t[0].start))
    assert len(clusters) >= 5  # centroid rule caps the diameter
    assert max(len(c.members) for c in clusters) < 15


def test_precise_needs_tight_breakpoints():
    cluster = SiteCluster(svtype="DEL", members=[_mc("A", 10_000, 60_000)])
    assert _compatible(_mc("B", 10_050, 60_100), cluster)          # within tolerance
    assert not _compatible(_mc("C", 14_000, 64_000), cluster)      # RO ok, bp too far


def test_hwe_extremes():
    # perfect equilibrium-ish: common variant
    assert hwe_exact(50, 25, 25) > 0.01
    # gross excess of hets (typical artifact signature)
    assert hwe_exact(100, 0, 0) < 1e-6
    # monomorphic
    assert hwe_exact(0, 0, 100) == 1.0


def _hwe_brute_force(het, hom_alt, hom_ref):
    """Exact conditional test by direct enumeration (Wigginton definition)."""
    from math import exp, lgamma

    n = het + hom_alt + hom_ref
    rare = 2 * min(hom_alt, hom_ref) + het

    def log_prob(h):
        if (rare - h) % 2:
            return None
        r = (rare - h) // 2
        c = n - h - r
        if r < 0 or c < 0:
            return None
        # P(h | n, rare) ∝ n! / (r! h! c!) * 2^h
        return (
            lgamma(n + 1) - lgamma(r + 1) - lgamma(h + 1) - lgamma(c + 1)
            + h * np.log(2)
        )

    logs = {h: lp for h in range(rare + 1) if (lp := log_prob(h)) is not None}
    mx = max(logs.values())
    probs = {h: exp(lp - mx) for h, lp in logs.items()}
    total = sum(probs.values())
    obs = probs[het]
    return sum(p for p in probs.values() if p <= obs * 1.0000001) / total


def test_hwe_matches_brute_force():
    for het, hom_alt, hom_ref in [(5, 2, 13), (57, 13, 298), (10, 0, 90), (3, 3, 3)]:
        assert hwe_exact(het, hom_alt, hom_ref) == pytest.approx(
            _hwe_brute_force(het, hom_alt, hom_ref), rel=1e-6
        )
