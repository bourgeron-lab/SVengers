"""Filter-stage joining of depth calls left in abutting pieces (genotype/bridge.py)."""

import numpy as np
import pysam

from svengers.config import Config
from svengers.extract.store import EvidenceStoreWriter
from svengers.genotype.bridge import bridge_calls, join_runs


def _store(path):
    """chr1, 3 Mb at 1 kb: a heterozygous deletion 1.0-2.0 Mb with a 2 kb
    homozygous hole at 1.5 Mb (a common deletion on the remaining allele)."""
    w = EvidenceStoreWriter(path, sample="S1", reference_id="x", bin_size=1000)
    w.set_meta(mean_depth=30.0, coverage=30.0, ploidy={"chr1": 2})
    depth = np.full(3000, 30.0, np.float32)
    depth[1000:2000] = 15.0
    depth[1500:1502] = 0.0
    w.write_depth("chr1", depth, depth, np.zeros(3000, bool))
    pos = np.arange(10_000, 3_000_000, 7_000, dtype=np.uint32)
    vaf = np.where((pos >= 1_000_000) & (pos < 2_000_000), 1.0, 0.5).astype(np.float32)
    w.write_baf("chr1", pos, vaf, np.full(len(pos), 30, np.uint16))


def _calls(path, records, legacy=False):
    """records: (start, end, svtype, source, cn, rd), internal 0-based half-open."""
    hdr = pysam.VariantHeader()
    if legacy:
        hdr.add_line("##source=vedette-1.1.0")       # pre-1.4 coordinates
    else:
        hdr.add_line("##svengers_coords=padded")
    hdr.contigs.add("chr1", length=3_000_000)
    for a in ("DEL", "DUP"):
        hdr.add_line(f'##ALT=<ID={a},Description="{a}">')
    for k, n, t in (("SVTYPE", "1", "String"), ("END", "1", "Integer"), ("SVLEN", "1", "Integer"),
                    ("SOURCE", "1", "String"), ("NBINS", "1", "Integer"), ("MASKFRAC", "1", "Float"),
                    ("HOMLEN", "1", "Integer"),
                    ("CIPOS", "2", "Integer"), ("CIEND", "2", "Integer")):
        hdr.info.add(k, n, t, k)
    hdr.info.add("IMPRECISE", "0", "Flag", "IMPRECISE")
    for f in ("NoEvidence", "MaxN", "LowPosterior", "MinSize"):
        hdr.filters.add(f, None, None, f)
    for k, t in (("GT", "String"), ("GQ", "Integer"), ("CN", "Integer"), ("CNQ", "Float"),
                 ("RD", "Float"), ("PE", "Integer"), ("SR", "Integer"), ("CG", "Integer"),
                 ("BAFN", "Integer"), ("BAFDEV", "Float")):
        hdr.formats.add(k, "1", t, k)
    hdr.add_sample("S1")
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for i, (s, e, sv, src, cn, rd, *flt) in enumerate(records):
            rec = vf.new_record(contig="chr1", start=s if legacy else s - 1,
                                stop=e + 1 if legacy else e, alleles=("N", f"<{sv}>"),
                                id=f"svengers_{sv}_{i}")
            rec.info["SVTYPE"], rec.info["SOURCE"] = sv, src
            rec.info["SVLEN"] = (e - s) * (1 if sv == "DUP" else -1)
            rec.info["NBINS"] = -(-(e - s) // 1000)
            rec.info["MASKFRAC"] = 0.0
            rec.info["CIPOS"], rec.info["CIEND"] = (-1000, 1000), (-1000 - i, 1000 + i)
            for f in (flt[0] if flt else ()):
                rec.filter.add(f)
            smp = rec.samples["S1"]
            smp["GT"], smp["GQ"], smp["CN"], smp["RD"] = (0, 1), 50, cn, rd
            smp["PE"] = smp["SR"] = smp["CG"] = 0 if src == "DEPTH" else 5
            vf.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _read(path):
    with pysam.VariantFile(str(path)) as vf:
        has = "BRIDGED" in vf.header.info
        return [(r.id, r.start, r.stop, r.info.get("BRIDGED") if has else None,
                 dict(r.samples[0]), {**dict(r.info), "FILTER": sorted(r.filter.keys())}) for r in vf]


def test_a_deletion_split_by_a_homozygous_hole_is_joined(tmp_path):
    """The two CN1 pieces around the 2 kb CN0 hole are one call again (as
    `call` has made it since v1.12.0); the hole's own junction call stays."""
    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vin, vout = tmp_path / "in.vcf.gz", tmp_path / "out.vcf.gz"
    _calls(vin, [(1_000_000, 1_500_000, "DEL", "DEPTH", 1, 0.5),
                 (1_499_700, 1_502_100, "DEL", "BOTH", 0, 0.0),
                 (1_502_000, 2_000_000, "DEL", "DEPTH", 1, 0.5)])
    c = bridge_calls(vin, store, vout, Config())
    assert (c["calls"], c["joined"], c["absorbed"]) == (3, 1, 1)
    recs = _read(vout)
    assert [r[0] for r in recs] == ["svengers_DEL_0", "svengers_DEL_1"]
    rid, start, stop, n, fmt, info = recs[0]
    assert (start + 1, stop, n) == (1_000_000, 2_000_000, 2)     # spec coordinates
    assert info["SVLEN"] == -1_000_000 and info["NBINS"] == 1000 and info["CIEND"] == (-1002, 1002)
    assert fmt["CN"] == 1 and tuple(fmt["GT"]) == (0, 1) and abs(fmt["RD"] - 0.495) < 0.01
    assert recs[1][3] is None                                   # the nested call unchanged


def test_legacy_coordinates_are_kept(tmp_path):
    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vin, vout = tmp_path / "in.vcf.gz", tmp_path / "out.vcf.gz"
    _calls(vin, [(1_000_000, 1_500_000, "DEL", "DEPTH", 1, 0.5),
                 (1_502_000, 2_000_000, "DEL", "DEPTH", 1, 0.5)], legacy=True)
    assert bridge_calls(vin, store, vout, Config())["joined"] == 1
    (_id, start, stop, n, _fmt, info), = _read(vout)
    assert (start, stop - 1, n, info["SVLEN"]) == (1_000_000, 2_000_000, 2, -1_000_000)


def test_junction_calls_and_distant_or_discordant_pieces_are_not_joined(tmp_path):
    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vin, vout = tmp_path / "in.vcf.gz", tmp_path / "out.vcf.gz"
    _calls(vin, [(1_000_000, 1_500_000, "DEL", "BOTH", 1, 0.5),     # junction-supported
                 (1_502_000, 1_700_000, "DEL", "DEPTH", 1, 0.5),
                 (1_710_000, 1_800_000, "DEL", "DEPTH", 1, 0.5),    # gap 10 bins
                 (1_802_000, 2_000_000, "DEL", "DEPTH", 1, 0.8)])   # ratio differs by 0.3
    c = bridge_calls(vin, store, vout, Config())
    assert (c["joined"], c["absorbed"]) == (0, 0)
    assert len(_read(vout)) == 4


def test_join_runs_chain_and_short_side_fraction():
    """Joins chain; the gap must also be <= 5 % of the shorter side."""
    pieces = [(0, 100_000, 1, 0.5), (102_000, 200_000, 1, 0.5), (201_000, 300_000, 1, 0.52),
              (302_000, 330_000, 1, 0.5),       # 2 bins > 5 % of 28 bins
              (100_500, 101_500, 0, 0.0)]       # other CN: ignored
    assert join_runs(pieces, 1000, 3, 0.05, 0.15) == [[0, 1, 2]]
    assert join_runs(pieces, 1000, 0, 0.05, 0.15) == []


def test_filter_scores_the_joined_call(tmp_path):
    """`svengers filter` joins before scoring and leaves no temporary file."""
    from svengers.filter.pipeline import run_filter

    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vin, vout = tmp_path / "S1.svengers.vcf.gz", tmp_path / "S1.filtered.vcf.gz"
    _calls(vin, [(1_000_000, 1_500_000, "DEL", "DEPTH", 1, 0.5),
                 (1_502_000, 2_000_000, "DEL", "DEPTH", 1, 0.5)])
    counts = run_filter(vin, store, vout, Config())      # no tracks, no model
    assert counts["total"] == 1
    (_id, start, stop, n, _fmt, _info), = _read(vout)
    assert (start + 1, stop, n) == (1_000_000, 2_000_000, 2)
    assert sorted(p.name for p in tmp_path.iterdir() if "bridged" in p.name) == []


def test_masked_pieces_are_not_joined_and_tags_are_the_joined_calls(tmp_path):
    """A MaxN piece (mostly masked: a segmental duplication) is never joined, so
    the event does not run into the repeat; the joined call gets its own
    call-time tags, not those of its first piece."""
    store = tmp_path / "S1.evidence.zarr"
    _store(store)
    vin, vout = tmp_path / "in.vcf.gz", tmp_path / "out.vcf.gz"
    _calls(vin, [(850_000, 998_000, "DEL", "DEPTH", 1, 0.5, ("MaxN", "LowPosterior")),
                 (1_000_000, 1_400_000, "DEL", "DEPTH", 1, 0.5, ("LowPosterior",)),
                 (1_402_000, 2_000_000, "DEL", "DEPTH", 1, 0.5)])
    assert bridge_calls(vin, store, vout, Config())["joined"] == 1
    recs = _read(vout)
    assert [(r[1] + 1, r[2], r[3]) for r in recs] == [(850_000, 998_000, None), (1_000_000, 2_000_000, 2)]
    assert recs[0][5]["FILTER"] == ["LowPosterior", "MaxN"]
    assert recs[1][5]["FILTER"] == ["PASS"]          # LowPosterior on 40 % of the bins: not inherited
