"""Sex-chromosome karyotype from the evidence store and the PoN profiles."""

import numpy as np
import pytest
import zarr

from vedette.extract.store import EvidenceStore, EvidenceStoreWriter, _create_array
from vedette.karyotype import classify, is_anomaly, sex_chromosome_stats

BIN = 10_000
N = 3_000   # bins per contig (30 Mb)


def _pon(tmp_path):
    """chrX: bins 0-199 pseudo-autosomal (1.0 in both sexes), the rest 1.0 in
    women / 0.5 in men; chrY: bins 0-999 a repeat seen in women too."""
    root = zarr.open(str(tmp_path / "pon.zarr"), mode="w")
    for c in ("chrX", "chrY"):
        g = root.create_group(c)
        if c == "chrX":
            f = np.ones(N)
            m = np.full(N, 0.5)
            m[:200] = 1.0
        else:
            f = np.zeros(N)
            f[:1000] = 0.4
            m = np.full(N, 0.5)
        for name, arr in (("median_female", f), ("median_male", m),
                          ("mad_female", np.full(N, 0.05)), ("mad_male", np.full(N, 0.05))):
            _create_array(g, name, arr, N)
    return root


def _store(tmp_path, x_copies, y_copies, het_x=None, x_baf_shift=0.0, seed=0):
    rng = np.random.default_rng(seed)
    path = tmp_path / f"S{x_copies}{y_copies}.evidence.zarr"
    w = EvidenceStoreWriter(path, sample="S", reference_id="x", bin_size=BIN)
    for c, cn in (("chr2", 2), ("chr7", 2), ("chr12", 2), ("chrX", x_copies), ("chrY", y_copies)):
        depth = (15.0 * cn + rng.normal(0, 1.0, N)).clip(0).astype(np.float32)
        if c == "chrX":
            depth[:200] = 30.0                      # pseudo-autosomal: 2 copies in both sexes
        w.write_depth(c, depth, depth, np.zeros(N, dtype=bool))
        n = 20_000
        pos = np.sort(rng.integers(0, N * BIN, n)).astype(np.uint32)
        het_rate = 0.6 if c != "chrX" else (het_x if het_x is not None else (0.6 if x_copies >= 2 else 0.003))
        het = rng.random(n) < het_rate
        side = rng.choice([-1, 1], n)
        shift = x_baf_shift if c == "chrX" else 0.0
        vaf = np.where(het, 0.5 + side * shift + rng.normal(0, 0.07, n), 1.0)
        w.write_baf(c, pos, np.clip(vaf, 0.0, 1.0).astype(np.float32), np.zeros(n, np.uint16))
    return EvidenceStore(path)


@pytest.mark.parametrize("x,y,expected", [(2, 0, "XX"), (1, 1, "XY"), (1, 0, "X0"),
                                          (2, 1, "XXY"), (1, 2, "XYY")])
def test_karyotypes(tmp_path, x, y, expected):
    st = sex_chromosome_stats(_store(tmp_path, x, y), _pon(tmp_path))
    assert st["x_cn"] == pytest.approx(x, abs=0.1)
    assert (st["y_cn"] or 0) == pytest.approx(y, abs=0.1)
    k = classify(st)
    assert k["karyotype"] == expected
    assert is_anomaly(k) == (expected not in ("XX", "XY"))


def test_triple_x_and_its_baf(tmp_path):
    st = sex_chromosome_stats(_store(tmp_path, 3, 0, x_baf_shift=1 / 6), _pon(tmp_path))
    assert classify(st)["karyotype"] == "XXX"
    assert st["x_baf"] == pytest.approx(1 / 6, abs=0.03)


def test_pedigree_discordance_contamination_and_mosaic():
    xy = {"x_cn": 1.0, "y_cn": 0.98, "x_het": 0.01, "x_baf": 0.38}
    assert "pedigree" in classify(xy, "F")["note"]
    assert classify(xy, "M") == {"karyotype": "XY", "note": ""}
    assert "contamination" in classify({**xy, "x_het": 0.4}, "M")["note"]
    assert classify({"x_cn": 1.7, "y_cn": 0.0, "x_het": 1.0, "x_baf": 0.06})["karyotype"] == "intermédiaire"
    assert "perte partielle du Y" in classify({**xy, "y_cn": 0.55})["note"]


@pytest.mark.parametrize("kar,expected", [
    ("XX", set()), ("XY", set()), ("?", set()), ("intermédiaire", set()),
    ("XXY", {"X"}), ("XXX", {"X"}), ("X0", {"X"}), ("XYY", {"Y"}), ("XXYY", {"X", "Y"}),
])
def test_aneuploid_sex_chromosomes(kar, expected):
    from vedette.karyotype import aneuploid_sex_chromosomes

    assert aneuploid_sex_chromosomes(kar) == expected
