import numpy as np

from svengers.breakpoints.cluster import cluster_cigar_dels
from svengers.breakpoints.resolve import merge_junction_sets
from svengers.config import BreakpointConfig


def _cig(entries):
    return {
        "pos": np.array([e[0] for e in entries], np.uint32),
        "svlen": np.array([e[1] for e in entries], np.uint32),
    }


def test_concordant_cigar_dels_form_junction():
    cig = _cig([(10_000, 120)] * 5 + [(10_002, 118)] * 2)
    junctions = cluster_cigar_dels(cig, BreakpointConfig())
    assert len(junctions) == 1
    j = junctions[0]
    assert j.svtype == "DEL" and j.precise
    assert j.cig == 7
    assert abs(j.start - 10_000) <= 2 and abs(j.size - 120) <= 2


def test_below_support_threshold_dropped():
    cig = _cig([(10_000, 120)] * 2)  # min_cigar_support = 3
    assert cluster_cigar_dels(cig, BreakpointConfig()) == []


def test_distinct_dels_not_merged():
    cig = _cig([(10_000, 120)] * 3 + [(10_000, 400)] * 3 + [(50_000, 120)] * 3)
    junctions = cluster_cigar_dels(cig, BreakpointConfig())
    assert len(junctions) == 3
    sizes = sorted(j.size for j in junctions)
    assert sizes[0] == 120 and sizes[-1] == 400


def test_small_dels_below_min_size_dropped():
    cig = _cig([(10_000, 40)] * 5)  # < min_size (50)
    assert cluster_cigar_dels(cig, BreakpointConfig()) == []


def test_dup_junction_only_requires_sr():
    from svengers.breakpoints.cluster import Junction
    from svengers.config import Config
    from svengers.merge_candidates import reconcile

    cfg = Config()
    pe_only_dup = Junction("DUP", 10_000, 10_500, (0, 0), (0, 0), 6, 0, 0, False)
    sr_dup = Junction("DUP", 30_000, 30_500, (0, 0), (0, 0), 4, 3, 3, True)
    pe_only_del = Junction("DEL", 50_000, 50_500, (0, 0), (0, 0), 6, 0, 0, False)
    calls = reconcile("chr1", [], [pe_only_dup, sr_dup, pe_only_del], 1000, 2, cfg)
    kinds = {(c.svtype, c.start) for c in calls}
    assert ("DUP", 30_000) in kinds          # SR-supported DUP kept
    assert ("DUP", 10_000) not in kinds      # everted-only DUP dropped
    assert ("DEL", 50_000) in kinds          # DEL rule unchanged


def test_cig_junction_merges_into_pe_junction():
    from svengers.breakpoints.cluster import Junction

    pe_junc = Junction("DEL", 9_990, 10_115, (-50, 50), (-50, 50), 4, 0, 0, False)
    cig_junc = cluster_cigar_dels(_cig([(10_000, 120)] * 5), BreakpointConfig())
    merged = merge_junction_sets([pe_junc], cig_junc)
    assert len(merged) == 1
    j = merged[0]
    assert j.precise and j.cig == 5 and j.pe == 4
    assert j.start == 10_000 and j.end == 10_120
