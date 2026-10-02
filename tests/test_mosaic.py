"""Mosaic chromosomal-scale CNVs from the evidence store."""

import numpy as np
import pytest

from vedette.extract.store import EvidenceStore, EvidenceStoreWriter
from vedette.mosaic import expected_baf_excess, read_tsv, scan_store, write_tsv

BIN = 10_000            # 100 bins per 1 Mb window: fast, same window logic
CONTIG_MB = 40


def _store(tmp_path, events, n_contigs=12, seed=0):
    """Synthetic store: `events` = {contig: (start_mb, end_mb, cn, baf_shift)}.

    baf_shift=True gives the segment the allelic imbalance its cell fraction
    implies; False keeps BAF at 0.5 (a depth artefact)."""
    rng = np.random.default_rng(seed)
    path = tmp_path / "S.evidence.zarr"
    w = EvidenceStoreWriter(path, sample="S", reference_id="x", bin_size=BIN)
    n_bins = CONTIG_MB * 1_000_000 // BIN
    for k in range(1, n_contigs + 1):
        c = f"chr{k}"
        cn = np.full(n_bins, 2.0)
        imb = np.zeros(n_bins)
        evs = events.get(c, [])
        for a, b, ev_cn, shift in ([evs] if isinstance(evs, tuple) else evs):
            sl = slice(a * 1_000_000 // BIN, b * 1_000_000 // BIN)
            cn[sl] = ev_cn
            if shift:
                f = min(abs(ev_cn - 2), 0.9)
                imb[sl] = expected_baf_excess("DEL" if ev_cn < 2 else "DUP", f)
        corr = (30.0 * cn / 2 + rng.normal(0, 1.5, n_bins)).astype(np.float32)
        w.write_depth(c, corr, corr, np.zeros(n_bins, dtype=bool))
        pos = np.sort(rng.integers(0, CONTIG_MB * 1_000_000, 40 * CONTIG_MB)).astype(np.uint32)
        side = rng.choice([-1, 1], len(pos))
        vaf = 0.5 + side * imb[pos // BIN] + rng.normal(0, 0.04, len(pos))
        w.write_baf(c, pos, np.clip(vaf, 0.12, 0.88).astype(np.float32),
                    np.zeros(len(pos), np.uint16))
    return EvidenceStore(path)


def test_mosaic_loss_and_gain_are_found_with_their_fraction(tmp_path):
    """A two-thirds clone: 10q-like loss at CN 1.33 and a whole-arm gain at
    CN 2.67, both with the allelic imbalance that fraction implies."""
    st = _store(tmp_path, {"chr3": (10, 40, 1.33, True), "chr5": (0, 25, 2.67, True)})
    calls = {c.contig: c for c in scan_store(st)}
    assert set(calls) == {"chr3", "chr5"}
    loss, gain = calls["chr3"], calls["chr5"]
    assert loss.svtype == "DEL" and gain.svtype == "DUP"
    assert loss.fraction == pytest.approx(0.67, abs=0.05)
    assert gain.fraction == pytest.approx(0.67, abs=0.05)
    assert loss.mosaic and gain.mosaic
    assert abs(loss.start - 10_000_000) <= 1_000_000 and loss.end == 40_000_000


def test_depth_shift_without_allelic_imbalance_is_rejected(tmp_path):
    """GC/wave artefacts and pericentromeric repeats move depth, not BAF."""
    st = _store(tmp_path, {"chr3": (10, 40, 1.33, False)})
    assert scan_store(st) == []


def test_normal_genome_is_silent_and_tsv_round_trips(tmp_path):
    st = _store(tmp_path, {})
    assert scan_store(st) == []
    (tmp_path / "b").mkdir()
    calls = scan_store(_store(tmp_path / "b", {"chr4": (5, 30, 1.4, True)}))
    write_tsv(calls, tmp_path / "m.tsv")
    rows = read_tsv(tmp_path / "m.tsv")
    assert len(rows) == 1 and rows[0]["contig"] == "chr4" and rows[0]["mosaic"] is True


def test_mosaic_flags_on_family_rows():
    from vedette.report import mosaic_flags

    seg = [{"contig": "chr10", "start": 83_000_000, "end": 133_000_000, "svtype": "DEL",
            "fraction": 0.67, "mosaic": True}]
    inside = {"contig": "chr10", "pos": 101_312_001, "end": 103_463_001, "svtype": "DEL"}
    other = {"contig": "chr2", "pos": 1_000, "end": 9_000, "svtype": "DEL"}
    assert mosaic_flags(seg, inside) == ["mosaïque 67 % chr10:83-133 Mb"]
    assert mosaic_flags(seg, other) == ["porteur mosaïque (1 segment)"]
    assert mosaic_flags([dict(seg[0], mosaic=False)], inside) == []
    assert mosaic_flags([], inside) == []


def test_blend_of_integer_levels_is_not_a_mosaic(tmp_path):
    """PMS: a constitutional 16p13 duplication beside a deleted window (2.9,
    2.5, 3.0, 1.15) and a 22q13 deletion with an intermediate edge (1.66,
    1.65, 1.03, 1.02) have intermediate medians but are not clones."""
    import numpy as np

    from vedette.mosaic import homogeneous

    assert not homogeneous(np.array([2.9, 2.5, 3.0, 1.15]), 2.5)
    assert not homogeneous(np.array([1.66, 1.65, 1.03, 1.02]), 1.65)
    assert homogeneous(np.array([1.33, 1.35, 1.33, 1.25, 1.33]), 1.33)
    # and a full-copy block resolved by the segmenter is a gain, not a mosaic
    st = _store(tmp_path, {"chr6": [(10, 13, 3.0, True), (13, 15, 1.0, True)]})
    chr6 = [c for c in scan_store(st) if c.contig == "chr6"]
    assert chr6 and not any(c.mosaic for c in chr6)


def test_short_terminal_mosaic_reaches_the_contig_end(tmp_path):
    """A PMS proband's terminal 22q13 loss at CN 1.62 over ~3 Mb, out of reach
    at a 5 Mb minimum; its end must stop at the contig, not at the window."""
    st = _store(tmp_path, {"chr8": (37, 40, 1.62, True)})
    calls = [c for c in scan_store(st) if c.contig == "chr8"]
    assert len(calls) == 1 and calls[0].mosaic
    assert calls[0].fraction == pytest.approx(0.38, abs=0.06)
    assert calls[0].end == CONTIG_MB * 1_000_000
