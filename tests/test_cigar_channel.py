import numpy as np

from vedette.breakpoints.cluster import cluster_cigar_dels
from vedette.breakpoints.resolve import merge_junction_sets
from vedette.config import BreakpointConfig


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


def test_cig_junction_merges_into_pe_junction():
    from vedette.breakpoints.cluster import Junction

    pe_junc = Junction("DEL", 9_990, 10_115, (-50, 50), (-50, 50), 4, 0, 0, False)
    cig_junc = cluster_cigar_dels(_cig([(10_000, 120)] * 5), BreakpointConfig())
    merged = merge_junction_sets([pe_junc], cig_junc)
    assert len(merged) == 1
    j = merged[0]
    assert j.precise and j.cig == 5 and j.pe == 4
    assert j.start == 10_000 and j.end == 10_120
