"""Inter-chromosomal junctions and terminal-loss mechanisms from the store."""

import numpy as np

from vedette.bnd import (
    BND_PROMISCUOUS_PARTNERS,
    BND_RECURRENT_MIN,
    Bnd,
    flag_junctions,
    loss_mechanism,
    scan_junctions,
    telomeric,
    terminal_losses,
)
from vedette.breakpoints.cluster import SRBreakpoint, cluster_sr
from vedette.config import BreakpointConfig
from vedette.extract.store import (
    PE_DEL,
    PE_INTER,
    PE_INV,
    SR_LEFT_CLIP,
    SR_RIGHT_CLIP,
    EvidenceStore,
    EvidenceStoreWriter,
)

NAMES = ["chr8", "chr22", "chrUn_GL000220v1"]
LENGTHS = {"chr8": 145_138_636, "chr22": 50_818_468, "chrUn_GL000220v1": 161_802}
TELO = b"TTAGGG" * 12


def _store(tmp_path, pairs=(), clips=None, name="S", inv=()):
    """`pairs`: (pos1 on chr8, pos2 on chr22, mapq) INTER pairs; `inv`:
    (pos1, pos2) same-strand pairs on chr22; `clips`: {contig: [(pos, side, seq)]}."""
    w = EvidenceStoreWriter(tmp_path / f"{name}.evidence.zarr", sample=name, reference_id="x",
                            bin_size=1000)
    w.set_meta(insert_mu=350.0, insert_sigma=50.0)
    pairs = list(pairs) + [(1_000_000 + 7919 * i, 2_000_000, 60) for i in range(5)]  # scattered
    p1 = np.array([p[0] for p in pairs])
    p2 = np.array([p[1] for p in pairs])
    w.write_pe("chr8", p1, p2, np.full(len(pairs), 1), np.full(len(pairs), PE_INTER),
               np.array([p[2] for p in pairs]))
    p1 = np.array([100] + [p[0] for p in inv])
    p2 = np.array([600] + [p[1] for p in inv])
    w.write_pe("chr22", p1, p2, np.full(len(p1), 21), np.array([PE_DEL] + [PE_INV] * len(inv)),
               np.full(len(p1), 60))
    for c in ("chr8", "chr22"):
        rows = (clips or {}).get(c, [])
        w.write_sr(c, np.array([r[0] for r in rows], dtype=np.int64),
                   np.array([r[1] for r in rows], dtype=np.int64),
                   np.full(len(rows), 60), np.array([len(r[2]) for r in rows], dtype=np.int64),
                   [r[2] for r in rows])
    return EvidenceStore(tmp_path / f"{name}.evidence.zarr")


def _translocation(n=12, mapq=60):
    """der(22)t(8;22): chr8 joined at 119,479,700, chr22 at 47,821,250."""
    return [(119_479_400 + 20 * i, 47_821_300 + 15 * i, mapq) for i in range(n)]


def test_junction_found_and_placed_on_its_split_reads(tmp_path):
    clips = {"chr8": [(119_479_750, SR_LEFT_CLIP, b"ACGT" * 8)] * 4,
             "chr22": [(47_821_250, SR_RIGHT_CLIP, b"TTGA" * 8)] * 3}
    st = _store(tmp_path, _translocation(), clips)
    bnds = scan_junctions(st, NAMES)
    assert len(bnds) == 1
    b = bnds[0]
    assert (b.contig_a, b.pos_a, b.side_a) == ("chr8", 119_479_750, "L")
    assert (b.contig_b, b.pos_b, b.side_b) == ("chr22", 47_821_250, "R")
    assert b.pe == 12 and b.sr_a == 4 and b.sr_b == 3


def test_weak_or_unmappable_pairs_make_no_junction(tmp_path):
    assert scan_junctions(_store(tmp_path, _translocation(n=3)), NAMES) == []
    (tmp_path / "q").mkdir()
    assert scan_junctions(_store(tmp_path / "q", _translocation(mapq=5)), NAMES) == []


def _bnd(s, ca, pa, cb, pb, pe=10):
    return Bnd(s, ca, pa, "?", cb, pb, "?", pe, 0, 0, 60.0)


def test_flags_promiscuous_recurrent_subtelomeric_unplaced():
    """A locus joined to many loci is a repeat (chr2:32.9 Mb with 22q13
    breakpoints all over 41-50 Mb); the same two loci joined in many genomes
    an artefact; a family's shared translocation is neither."""
    by = {f"S{i}": [_bnd(f"S{i}", "chr8", 32_900_000 + i, "chr22", 41_000_000 + 900_000 * i)]
          for i in range(BND_PROMISCUOUS_PARTNERS)}
    for i in range(BND_RECURRENT_MIN):
        by[f"S{i}"].append(_bnd(f"S{i}", "chr8", 4_795_714, "chr22", 48_985_660 + i))
    by["T"] = [_bnd("T", "chr8", 119_479_750, "chr22", 47_821_250),
               _bnd("T", "chr8", 145_000_000, "chr22", 46_000_000),
               _bnd("T", "chr22", 45_000_000, "chrUn_GL000220v1", 100_000)]
    by["T_parent"] = [_bnd("T_parent", "chr8", 119_479_700, "chr22", 47_821_300)]
    flag_junctions(by, LENGTHS)
    assert {b.flag for s in by if s.startswith("S") for b in by[s]} == {"promiscuous", "recurrent"}
    assert [b.flag for b in by["T"]] == ["", "subtelomeric", "unplaced"]
    assert by["T_parent"][0].flag == ""


def test_recurrent_same_strand_junction_with_support_is_a_common_inversion():
    """In many genomes, a same-strand junction with 12 pairs or more is a
    common inversion (17q21.31 H2); weaker ones and inter-chromosomal ones
    stay recurrent (artefacts)."""
    from vedette.vcfout import bnd_filter

    by = {}
    for i in range(BND_RECURRENT_MIN):
        inv = Bnd(f"S{i}", "chr17", 45_585_000, "R", "chr17", 46_265_000, "R", 30 if i else 5, 4, 4,
                  60.0, kind="INV")
        tra = Bnd(f"S{i}", "chr8", 4_795_714, "R", "chr22", 48_985_660, "L", 20, 4, 4, 60.0)
        by[f"S{i}"] = [inv, tra]
    flag_junctions(by, LENGTHS)
    assert [by[f"S{i}"][0].flag for i in range(2)] == ["recurrent", "common"]
    assert {by[s][1].flag for s in by} == {"recurrent"}
    assert bnd_filter(by["S1"][0]) == "Common"


def test_terminal_losses_reach_the_end_and_are_long():
    calls = [("chr22", 47_821_250, 50_818_000), ("chr22", 50_600_000, 50_602_000),
             ("chr22", 30_000_000, 31_000_000)]
    mosaic = [{"contig": "chr22", "start": 41_993_000, "end": 50_819_000, "svtype": "DEL"}]
    assert terminal_losses(calls, mosaic, LENGTHS) == [
        ("chr22", 47_821_250, 50_818_000, "call"), ("chr22", 41_993_000, 50_819_000, "mosaic")]


def test_a_loss_called_in_pieces_breaks_at_the_inner_one():
    calls = [("chr22", 43_570_000, 49_320_000), ("chr22", 49_322_000, 50_818_468)]
    pieces = [("chr22", 41_000_000, 41_100_000), ("chr22", 42_000_000, 42_975_000),  # LowSQ
              ("chr22", 42_980_000, 42_995_000), ("chr22", 43_000_000, 43_560_000)]
    assert terminal_losses(calls, [], LENGTHS, pieces=pieces) == [
        ("chr22", 43_000_000, 50_818_468, "call")]


def test_a_loss_with_callable_sequence_past_it_is_not_terminal():
    """The common 2q37.3 deletion ends ~100 kb before chr2's end, where the
    reference is still callable; 22q's last 19 kb are masked."""
    lengths = {"chr2": 242_193_529, "chr22": 50_818_468}
    callable_at = {"chr2": np.arange(0, 242_159_000, 1000),
                   "chr22": np.arange(0, 50_800_000, 1000)}.get
    calls = [("chr2", 241_911_101, 242_093_628), ("chr22", 49_322_000, 50_799_001)]
    assert terminal_losses(calls, [], lengths, callable_at) == [
        ("chr22", 49_322_000, 50_799_001, "call")]


def test_mechanisms_translocation_telomere_unknown(tmp_path):
    loss = ("chr22", 47_821_250, 50_818_000, "call")
    trans = [_bnd("S", "chr8", 119_479_750, "chr22", 47_821_250, pe=16)]
    st = _store(tmp_path, clips={"chr22": [(47_821_250, SR_RIGHT_CLIP, TELO)] * 6})
    m = loss_mechanism(st, "S", loss, trans, [("chr8", 119_477_000, 145_000_000)], LENGTHS)
    assert m.mechanism == "translocation" and m.partner == "chr8:119479750"
    assert m.partner_gain == "chr8:119477000-145000000" and m.junction_pe == 16
    # no junction: the telomere repeats clipped on the retained side
    m = loss_mechanism(st, "S", loss, [], [], LENGTHS)
    assert m.mechanism == "telomere" and m.telo_clips == 6
    # telomere repeats clipped on the lost side are not a healed telomere
    (tmp_path / "b").mkdir()
    st = _store(tmp_path / "b", clips={"chr22": [(47_821_250, SR_LEFT_CLIP, TELO)] * 6})
    assert loss_mechanism(st, "S", loss, [], [], LENGTHS).mechanism == "unknown"
    # a junction to a repeat is no translocation: at the breakpoint itself it
    # says the end went onto unplaced sequence (rDNA, acrocentric short arms)
    junk = Bnd(**dict(vars(trans[0]), flag="promiscuous"))
    m = loss_mechanism(st, "S", loss, [junk], [], LENGTHS)
    assert m.mechanism == "repeat" and m.partner == "chr8:119479750"
    far = Bnd(**dict(vars(junk), pos_b=47_826_250))     # 5 kb off
    assert loss_mechanism(st, "S", loss, [far], [], LENGTHS).mechanism == "unknown"


def test_telomeric_clip():
    assert telomeric(b"GGTG" + TELO) and telomeric(b"CCCTAA" * 4)
    assert not telomeric(b"TTAGGGTTAGGGACGT")


def _fold(n=6, lo=49_329_800, hi=49_347_900):
    """++ pairs only: a fold-back joining 49,330,000 to 49,348,000 inverted."""
    return [(lo + 30 * i, hi + 25 * i) for i in range(n)]


def test_inversion_junctions_reciprocal_or_fold_back(tmp_path):
    clips = {"chr22": [(49_330_000, SR_RIGHT_CLIP, b"ACGT" * 8)] * 4}
    st = _store(tmp_path, clips=clips, inv=_fold() + [(49_300_000, 49_300_400)] * 5)  # + short
    (b,) = scan_junctions(st, NAMES)
    assert (b.kind, b.contig_a, b.pos_a, b.contig_b, b.pe) == ("FOLD", "chr22", 49_330_000, "chr22", 6)
    # -- pairs right of the first breakpoint too: a reciprocal inversion
    (tmp_path / "r").mkdir()
    st = _store(tmp_path / "r", clips=clips, inv=_fold() + _fold(4, 49_330_100, 49_348_100))
    (b,) = scan_junctions(st, NAMES)
    assert b.kind == "INV" and b.pe == 10


def test_inv_dup_del_is_a_fold_back_next_to_the_loss():
    """der(22): 41.5-49.33 Mb duplicated inverted, folded back at 49.33-49.35,
    the rest of 22q lost."""
    loss = ("chr22", 49_348_000, 50_818_468, "call")
    fold = Bnd("S", "chr22", 49_330_000, "R", "chr22", 49_348_050, "R", 5, 3, 2, 60.0, kind="FOLD")

    class _NoClips:
        def sr(self, *a, **k):
            return {"pos": np.array([], dtype=np.int64), "side": np.array([], dtype=np.int64),
                    "seq": []}

    m = loss_mechanism(_NoClips(), "S", loss, [fold], [("chr22", 41_536_000, 49_330_000)], LENGTHS)
    assert m.mechanism == "inv-dup" and m.partner == "chr22:49330000"
    assert m.partner_gain == "chr22:41536000-49330000" and m.junction_pe == 5
    # a short spacer: both ends within the window, the one at the breakpoint is its own
    short = Bnd("S", "chr22", 47_414_330, "L", "chr22", 47_420_204, "R", 19, 8, 9, 60.0, kind="FOLD")
    m = loss_mechanism(_NoClips(), "S", ("chr22", 47_420_000, 50_818_468, "call"), [short], [],
                       LENGTHS)
    assert m.mechanism == "inv-dup" and m.partner == "chr22:47414330"


def test_a_parent_with_the_junction_carries_the_balanced_translocation():
    from vedette.bnd import junction_carriers, terminal_label

    row = {"contig": "chr22", "start": 47_821_250, "end": 50_818_468, "partner": "chr8:119479750",
           "partner_gain": "chr8:119477000-145138636", "mechanism": "translocation"}
    jn = {"contig_a": "chr8", "pos_a": "119479700", "contig_b": "chr22", "pos_b": "47821300"}
    other = {"contig_a": "chr8", "pos_a": "119479700", "contig_b": "chr22", "pos_b": "30000000"}
    found = junction_carriers(row, {"father": [other], "mother": [other, jn]}, LENGTHS)
    assert found == ["mother"]
    assert terminal_label(row, found) == (
        "translocation chr8:119,479,750, gain chr8:119.5-145.1 Mb, jonction chez mother")
    assert junction_carriers(dict(row, partner=""), {"mother": [jn]}, LENGTHS) == []


def test_a_parent_with_the_same_loss_transmitted_it():
    from vedette.bnd import loss_carriers, terminal_label

    row = {"contig": "chr5", "start": 0, "end": 141_001, "partner": "", "partner_gain": "",
           "mechanism": "telomere"}
    mother = [dict(row, start=0, end=141_001)]
    father = [dict(row, contig="chr22", start=49_000_000, end=50_818_468)]
    shorter = [dict(row, end=60_000)]                       # overlap below one half
    assert loss_carriers(row, {"father": father, "mother": mother, "other": shorter}) == ["mother"]
    assert terminal_label(row, [], ["mother"]) == "télomère néoformé, héritée de mother"


def test_terminal_flag_on_the_family_row_of_the_loss():
    from vedette.report import terminal_flags

    losses = [{"contig": "chr22", "start": 43_000_000, "end": 50_818_468,
               "label": "télomère néoformé"}]
    assert terminal_flags(losses, {"contig": "chr22", "pos": 49_322_001, "end": 50_818_468,
                                   "svtype": "DEL"}) == ["perte terminale : télomère néoformé"]
    assert terminal_flags(losses, {"contig": "chr22", "pos": 30_000_001, "end": 30_100_000,
                                   "svtype": "DEL"}) == []
    assert terminal_flags(losses, {"contig": "chr22", "pos": 45_000_001, "end": 45_100_000,
                                   "svtype": "DUP"}) == []


def _cluster_sr_loop(sr, config):
    """The per-clip loop cluster_sr replaced: the reference it must match."""
    out = []
    for side in (SR_LEFT_CLIP, SR_RIGHT_CLIP):
        sel = sr["side"] == side
        if not sel.any():
            continue
        pos = sr["pos"][sel].astype(np.int64)
        mapq = sr["mapq"][sel].astype(np.float64)
        clip = sr["cliplen"][sel].astype(np.int64)
        order = np.argsort(pos, kind="stable")
        pos, mapq, clip = pos[order], mapq[order], clip[order]
        start = 0
        for i in range(1, len(pos) + 1):
            if i == len(pos) or pos[i] - pos[i - 1] > config.sr_pos_tolerance:
                n = i - start
                if n >= config.min_sr:
                    out.append(SRBreakpoint(int(np.median(pos[start:i])), side, n,
                                            float(mapq[start:i].mean()), int(clip[start:i].max())))
                start = i
    out.sort(key=lambda b: b.pos)
    return out


def test_vectorised_cluster_sr_matches_the_loop():
    rng = np.random.default_rng(7)
    for _ in range(20):
        k = int(rng.integers(1, 3000))
        sr = {"pos": np.sort(rng.integers(0, 40_000, k)), "side": rng.integers(0, 2, k),
              "mapq": rng.integers(0, 61, k), "cliplen": rng.integers(5, 150, k)}
        cfg = BreakpointConfig(min_sr=int(rng.integers(1, 4)))
        assert cluster_sr(sr, cfg) == _cluster_sr_loop(sr, cfg)


def test_breakend_vcf(tmp_path):
    """der(22)t(8;22): chr22 kept left of its breakend (reads clipped on the
    right) joined to chr8 right of its own (clipped on the left)."""
    import pysam

    from vedette.reference import Contig
    from vedette.vcfout import write_bnd_vcf

    fa = tmp_path / "ref.fa"
    fa.write_text(">chr8\n" + "ACGT" * 500 + "\n>chr22\n" + "TTGCA" * 400 + "\n")
    pysam.faidx(str(fa))
    contigs = [Contig("chr8", 2000, "autosome"), Contig("chr22", 2000, "autosome")]
    bnds = [Bnd("S", "chr8", 1_750, "L", "chr22", 1_250, "R", 16, 5, 6, 57.9),
            Bnd("S", "chr8", 100, "?", "chr22", 200, "R", 4, 0, 3, 60.0),       # no orientation
            Bnd("S", "chr22", 300, "R", "chr22", 1_000, "R", 6, 3, 2, 60.0, kind="FOLD",
                flag="subtelomeric")]
    out = tmp_path / "S.bnd.vcf.gz"
    assert write_bnd_vcf(bnds, contigs, "S", out, fa) == 2
    with pysam.VariantFile(str(out)) as v:
        recs = {r.id: r for r in v}
    a, b = recs["vedette_BND_1_1"], recs["vedette_BND_1_2"]
    # chr8 left clip at 1,750 (0-based first aligned base): breakend base 1,751
    assert (a.chrom, a.pos, a.ref) == ("chr8", 1_751, "ACGT"[1_750 % 4])
    assert a.alts == (f"]chr22:1250]{a.ref}",) and a.info["MATEID"] == ("vedette_BND_1_2",)
    # chr22 right clip at 1,250 (first clipped base): last kept base 1,250
    assert (b.chrom, b.pos, b.ref) == ("chr22", 1_250, "TTGCA"[1_249 % 5])
    assert b.alts == (f"{b.ref}[chr8:1751[",) and b.samples["S"]["SR"] == 6
    assert list(a.filter) == ["PASS"]   # 16 pairs
    fold = recs["vedette_BND_2_1"]
    assert fold.alts == (f"{fold.ref}]chr22:1000]",) and fold.info["KIND"] == "FOLD"
    assert list(fold.filter) == ["Subtelomeric"]


def test_same_strand_junction_pass_unless_the_genome_is_noisy():
    from vedette.vcfout import BND_NOISY_SAME_STRAND, bnd_filter, noisy_genome

    inv = Bnd("S", "chr5", 1_000_000, "R", "chr5", 1_200_000, "R", 14, 4, 3, 60.0, kind="INV")
    weak = Bnd("S", "chr5", 3_000_000, "L", "chr5", 3_010_000, "L", 6, 2, 2, 60.0, kind="FOLD")
    assert not noisy_genome([inv, weak])
    assert bnd_filter(inv) == "PASS" and bnd_filter(weak) == "LowSupport"
    chimeric = [Bnd("S", "chr1", i * 10_000, "R", "chr1", i * 10_000 + 5_000, "R", 4, 1, 1, 60.0,
                    kind="FOLD") for i in range(BND_NOISY_SAME_STRAND + 1)]
    assert noisy_genome(chimeric + [inv])
    assert bnd_filter(inv, noisy=True) == "NoisyGenome"
    tra = Bnd("S", "chr5", 1_000_000, "R", "chr9", 1_200_000, "L", 14, 4, 3, 60.0)
    assert bnd_filter(tra, noisy=True) == "PASS"
