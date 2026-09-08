import numpy as np
import pytest

from vedette.breakpoints.cluster import Junction
from vedette.breakpoints.resolve import (
    _fit_align,
    merge_junction_sets,
    microhomology,
    recover_junctions,
)
from vedette.config import BreakpointConfig
from vedette.extract.store import SR_LEFT_CLIP, SR_RIGHT_CLIP


@pytest.fixture(scope="module")
def ref():
    rng = np.random.default_rng(123)
    return rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=10_000)


def _sr(entries):
    """entries: (pos, side, seq_bytes) — duplicated clips share a pileup."""
    pos = np.array([e[0] for e in entries], np.uint32)
    side = np.array([e[1] for e in entries], np.uint8)
    return (
        {
            "pos": pos,
            "side": side,
            "mapq": np.full(len(entries), 60, np.uint8),
            "cliplen": np.array([len(e[2]) for e in entries], np.uint16),
        },
        [e[2] for e in entries],
    )


def test_fit_align_exact_and_mismatch(ref):
    query = ref[3000:3060].copy()
    score, rs, re_ = _fit_align(query, ref[2000:4000])
    assert (rs, re_) == (1000, 1060)
    assert score == 2 * 60
    query[30] = ord("A") if query[30] != ord("A") else ord("C")
    score2, rs2, re2 = _fit_align(query, ref[2000:4000])
    assert (rs2, re2) == (1000, 1060)
    assert score2 == 2 * 59 - 3


def test_right_clip_downstream_is_del(ref):
    # DEL(2000, 2500): reads clipped at 2000, clip = sequence starting at 2500
    clip = bytes(ref[2500:2545])
    sr, seqs = _sr([(2000, SR_RIGHT_CLIP, clip), (2000, SR_RIGHT_CLIP, clip)])
    junctions = recover_junctions(sr, seqs, ref, BreakpointConfig())
    assert len(junctions) == 1
    j = junctions[0]
    assert (j.svtype, j.start, j.end) == ("DEL", 2000, 2500)
    assert j.precise and j.sr >= 4  # support counted on both breakpoint sides


def test_right_clip_upstream_is_tandem_dup(ref):
    # tandem DUP(4000, 4800): read crossing copy1|copy2 junction is clipped at
    # 4800 (dup end) and its clip is the dup START sequence
    clip = bytes(ref[4000:4042])
    sr, seqs = _sr([(4800, SR_RIGHT_CLIP, clip), (4800, SR_RIGHT_CLIP, clip)])
    junctions = recover_junctions(sr, seqs, ref, BreakpointConfig())
    assert len(junctions) == 1
    j = junctions[0]
    assert (j.svtype, j.start, j.end) == ("DUP", 4000, 4800)


def test_left_clip_geometry(ref):
    # DEL seen from the right anchor: LEFT clip at 2500, clip ENDS at 2000
    del_clip = bytes(ref[1955:2000])
    # DUP seen from the left anchor: LEFT clip at 4000, clip ENDS at 4800
    dup_clip = bytes(ref[4755:4800])
    sr, seqs = _sr([
        (2500, SR_LEFT_CLIP, del_clip), (2500, SR_LEFT_CLIP, del_clip),
        (4000, SR_LEFT_CLIP, dup_clip), (4000, SR_LEFT_CLIP, dup_clip),
    ])
    junctions = recover_junctions(sr, seqs, ref, BreakpointConfig())
    kinds = {(j.svtype, j.start, j.end) for j in junctions}
    assert ("DEL", 2000, 2500) in kinds
    assert ("DUP", 4000, 4800) in kinds


def test_noise_clip_no_junction(ref):
    rng = np.random.default_rng(7)
    noise = bytes(rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=45))
    sr, seqs = _sr([(3000, SR_RIGHT_CLIP, noise), (3000, SR_RIGHT_CLIP, noise)])
    assert recover_junctions(sr, seqs, ref, BreakpointConfig()) == []


def test_single_clip_below_support_threshold(ref):
    clip = bytes(ref[2500:2545])
    sr, seqs = _sr([(2000, SR_RIGHT_CLIP, clip), (2000, SR_RIGHT_CLIP, b"")])
    # only one usable clip (min support 2, and min_sr pileup needs the pair)
    assert recover_junctions(sr, seqs, ref, BreakpointConfig()) == []


def test_microhomology():
    base = np.frombuffer(b"AAAATTTTCCCCGGGGACGTACGTAAAATTTT", dtype=np.uint8)
    # ref[0:4] == ref[24:28] == AAAA and ref[4:8] == ref[28:32] == TTTT
    assert microhomology(base, 4, 28) == 8  # 4 forward (TTTT) + 4 backward (AAAA)
    rng = np.random.default_rng(1)
    rand = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=1000)
    assert microhomology(rand, 100, 600) <= 4  # random: tiny by chance


def test_merge_junction_sets_refines_and_appends():
    pe_junc = Junction(
        svtype="DUP", start=3990, end=4810, cipos=(-200, 200), ciend=(-200, 200),
        pe=6, sr_start=0, sr_end=0, precise=False,
    )
    resolved = [
        Junction("DUP", 4000, 4800, (0, 0), (0, 0), 0, 3, 3, True, homlen=7),
        Junction("DEL", 8000, 8300, (0, 0), (0, 0), 0, 2, 2, True, homlen=2),
    ]
    merged = merge_junction_sets([pe_junc], resolved)
    assert len(merged) == 2
    dup = next(j for j in merged if j.svtype == "DUP")
    assert (dup.start, dup.end, dup.precise, dup.homlen) == (4000, 4800, True, 7)
    assert dup.pe == 6  # PE support preserved
    assert any(j.svtype == "DEL" and j.start == 8000 for j in merged)
