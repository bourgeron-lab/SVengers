"""Per-call CNN tensors read from the evidence store."""

import numpy as np
import pysam
import pytest
import zarr

from vedette.extract.store import PE_DEL, SR_LEFT_CLIP, SR_RIGHT_CLIP, EvidenceStoreWriter, _create_array
from vedette.filter.tensor import CHANNELS, L, SCALARS, run_tensors

BIN = 1000
N = 300  # bins per contig (300 kb)
DEL = (100_000, 120_000)   # one-copy loss on chr1
C = {name: i for i, name in enumerate(CHANNELS)}


def _store(tmp_path, pon=None):
    rng = np.random.default_rng(0)
    w = EvidenceStoreWriter(tmp_path / "S.evidence.zarr", sample="S", reference_id="x", bin_size=BIN)
    w.set_meta(mean_depth=30.0, coverage=30.0, insert_mu=350.0, insert_sigma=70.0, insert_n_pairs=1000,
               ploidy={"chr1": 2, "chrX": 1}, wave_score=0.1, pon=str(pon) if pon else None)
    for c, base in (("chr1", 30.0), ("chrX", 15.0)):
        depth = (base + rng.normal(0, 0.5, N)).astype(np.float32)
        if c == "chr1":
            depth[DEL[0] // BIN: DEL[1] // BIN] /= 2
        w.write_depth(c, depth, depth, np.zeros(N, bool))
        pos = np.arange(500, N * BIN, 2_000)
        het = (pos < DEL[0]) | (pos >= DEL[1]) if c == "chr1" else np.zeros(len(pos), bool)
        w.write_baf(c, pos, np.where(het, 0.5, 1.0).astype(np.float32), np.full(len(pos), 30))
        # ten DEL-like pairs across the deletion, clips at both junctions
        n = 10 if c == "chr1" else 0
        w.write_pe(c, np.full(n, DEL[0] - 200), np.full(n, DEL[1] + 200), np.zeros(n), np.full(n, PE_DEL),
                   np.full(n, 60))
        w.write_sr(c, np.array([DEL[0]] * n + [DEL[1]] * n), np.array([SR_RIGHT_CLIP] * n + [SR_LEFT_CLIP] * n),
                   np.full(2 * n, 60), np.full(2 * n, 30))
        w.write_ins(c, np.full(N * 10, 350.0 * 20), np.full(N * 10, 20))
        w.write_cig(c, np.array([], np.uint32), np.array([], np.uint32))
    return tmp_path / "S.evidence.zarr"


def _pon(tmp_path):
    root = zarr.open(str(tmp_path / "pon.zarr"), mode="w")
    root.attrs["bin_size"] = BIN
    g = root.create_group("chr1")
    med = np.ones(N, np.float32)
    med[200:210] = 0.5          # a dip the whole panel shares
    _create_array(g, "median_all", med, N)
    _create_array(g, "mad_all", np.full(N, 0.05, np.float32), N)
    return tmp_path / "pon.zarr"


def _vcf(tmp_path, calls):
    h = pysam.VariantHeader()
    for c in ("chr1", "chrX"):
        h.contigs.add(c, length=N * BIN)
    for line in ('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">',
                 '##INFO=<ID=SOURCE,Number=1,Type=String,Description="s">',
                 '##INFO=<ID=SEGDUP_OVL,Number=1,Type=Float,Description="s">',
                 '##FILTER=<ID=LowSQ,Description="l">', '##FILTER=<ID=MaxN,Description="m">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">',
                 '##FORMAT=<ID=SQ,Number=1,Type=Float,Description="q">'):
        h.add_line(line)
    h.add_sample("S")
    path = tmp_path / "calls.vcf.gz"
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        for i, (contig, s, e, t, flt) in enumerate(calls):
            r = out.new_record(contig=contig, start=s, stop=e, alleles=("N", f"<{t}>"), id=f"c{i}")
            r.info["SVTYPE"], r.info["SOURCE"], r.info["SEGDUP_OVL"] = t, "DEPTH", 0.2
            r.filter.add(flt)
            r.samples["S"]["GT"], r.samples["S"]["SQ"] = (0, 1), 7.5
            out.write(r)
    return path


def _run(tmp_path, calls, pon=None, **kw):
    out = tmp_path / "t.npz"
    n = run_tensors(_vcf(tmp_path, calls), _store(tmp_path, pon), out, **kw)
    return n, np.load(out)


def test_deletion_shape(tmp_path):
    n, z = _run(tmp_path, [("chr1", *DEL, "DEL", "PASS"), ("chr1", 200_000, 200_500, "DEL", "PASS"),
                           ("chr1", 250_000, 270_000, "DEL", "MaxN")])
    assert n == 1                              # < 1 kb and hard-filtered calls dropped
    x = z["x"][0].astype(np.float32)
    assert x.shape == (len(CHANNELS), L) and z["s"].shape == (1, len(SCALARS))
    inside, flanks = x[C["in_call"]] > 0.99, x[C["in_call"]] == 0
    assert inside.sum() == pytest.approx(L / 3, abs=2)       # 20 kb call, 20 kb flanks
    assert np.median(x[C["depth"]][inside]) == pytest.approx(-1.0, abs=0.05)
    assert np.median(x[C["depth"]][flanks]) == pytest.approx(0.0, abs=0.05)
    assert x[C["het_density"]][inside].max() == 0 and x[C["het_density"]][flanks].max() > 0
    # pairs and clips sit at the two junctions, nowhere else
    for name in ("pe_del", "sr_left", "sr_right"):
        hits = np.nonzero(x[C[name]])[0]
        assert len(hits) and all(abs(i - L / 3) < 3 or abs(i - 2 * L / 3) < 3 for i in hits)
    assert np.allclose(x[C["ins_z"]], 0) and x[C["depth_masked"]].max() == 0
    assert not z["has_pon"] and np.allclose(x[C["pon_median"]], 0)
    assert z["call_id"][0] == "c0" and z["sq"][0] == pytest.approx(7.5)
    s = dict(zip(SCALARS, z["s"][0]))
    assert s["src_depth"] == 1 and s["segdup_ovl"] == pytest.approx(0.2) and s["ploidy"] == 1


def test_cohort_context_and_male_x(tmp_path):
    n, z = _run(tmp_path, [("chr1", 200_000, 210_000, "DEL", "LowSQ"),
                           ("chrX", 100_000, 120_000, "DUP", "LowSQ")], pon=_pon(tmp_path))
    assert n == 2 and z["has_pon"]
    x = z["x"].astype(np.float32)
    inside = x[0, C["in_call"]] > 0.99
    assert np.median(x[0, C["pon_median"]][inside]) == pytest.approx(-1.0, abs=0.01)  # panel dip
    assert np.median(x[0, C["pon_median"]][~inside]) == pytest.approx(0.0, abs=0.01)
    # a man's chrX at half the autosomal depth is his normal level
    assert np.median(x[1, C["depth"]]) == pytest.approx(0.0, abs=0.05)
    assert np.allclose(x[1, C["pon_median"]], 0)            # no chrX profile in this panel


def test_window_off_the_contig_start(tmp_path):
    _, z = _run(tmp_path, [("chr1", 2_000, 4_000, "DEL", "PASS")])
    x = z["x"][0].astype(np.float32)
    first = x[C["depth_masked"]][: L // 4]
    assert first.max() == 1.0 and x[C["depth_masked"]][-1] == 0.0
