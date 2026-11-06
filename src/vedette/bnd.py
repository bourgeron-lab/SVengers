"""Inter-chromosomal junctions and the mechanism of terminal losses.

A terminal deletion has one breakpoint, and what lies past it tells how the
chromosome was rescued: a new telomere added at the break (telomere healing),
or the end of another chromosome (an unbalanced translocation, whose partner
segment is then gained). This pass reads both from the evidence store alone:

1. inter-chromosomal pairs (PE_INTER, recorded once, in the store of the
   contig with the lower reference index) are clustered per contig pair, and
   each end is placed on its strongest split-read pileup when there is one;
   same-strand pairs (PE_INV) of 1 kb or more likewise make inversion
   junctions, reciprocal (pairs on both sides of the breakpoints) or
   fold-back (one side: one junction of an inverted duplication);
2. over the genomes of a run, a 10 kb locus joined to many distinct loci
   (BND_PROMISCUOUS_PARTNERS) is a repeat, and two loci joined in many
   genomes (BND_RECURRENT_*) an artefact or a common variant, not a
   breakpoint; a junction whose end is within SUBTEL_BP of a contig end lands
   in subtelomeric repeats, which is what a healed telomere's reads also hit;
3. each loss reaching a contig end (a PASS deletion or a mosaic segment)
   takes its mechanism from its breakpoint: a junction to another contig
   (translocation, unbalanced when the sample gains the partner's end), a
   fold-back on the retained side (inv-dup: the retained end duplicated
   inverted, folded back over a spacer), clipped reads made of TTAGGG
   repeats on the retained side or a subtelomeric junction (healed
   telomere), a junction to a repeat within 2 kb of the breakpoint (rDNA,
   acrocentric short arms: sequence the reference does not place),
   otherwise unknown. A loss called in abutting pieces breaks at
   the inner piece.

On the 63 PMS terminal 22q13 losses: ~20 healed telomeres (4-15 telomeric
clips each), translocations to 2q, 8q, 10p, 10q, 14q with 5-22 MAPQ-60 pairs
(two partner gains already seen by the mosaic scan sat exactly at the
junction); 10 losses seemed to break at 49.32 Mb, the outer piece of
losses called in two or more pieces from 43.6-49.1 Mb.
Pairs are recorded without read strands, so a junction's orientation comes
from the split reads; a junction without any has its ends on the pair
clusters, within an insert size.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path

import numpy as np

from .breakpoints.cluster import cluster_pe, cluster_sr
from .config import BreakpointConfig
from .extract.store import PE_INTER, PE_INV, SR_LEFT_CLIP, SR_RIGHT_CLIP

BND_MIN_PE = 4             # pairs to call a junction
BND_MIN_MAPQ = 20.0        # their mean MAPQ
LOCUS_BP = 10_000          # locus size for the run-wide counts
# a locus joined to this many distinct loci over the run is a repeat: chr2:32.9
# Mb pairs with ~24,000 loci in PMS (22q13 breakpoints at 41-50 Mb among them),
# a translocation's ends with 1-3. Counting genomes per locus does not work:
# PMS has ~12,700 raw junctions per genome, a median locus has some in 7.
BND_PROMISCUOUS_PARTNERS = 10
# the same two loci joined in this many genomes (or this share of the run) is
# a reference artefact or a common variant; a translocation is private or
# shared by a family (a balanced carrier and the unbalanced child)
BND_RECURRENT_MIN = 5
BND_RECURRENT_FRAC = 0.01
SUBTEL_BP = 500_000        # a junction end this close to a contig end is subtelomeric
TERMINAL_BP = 300_000      # a loss reaching this close to a contig end is terminal
TERMINAL_TAIL_BP = 10_000  # ... with at most this much callable sequence past it
CHAIN_GAP_BP = 20_000      # deletion pieces this close are one loss ...
CHAIN_MIN_BP = 50_000      # ... when this long (small noisy calls do not chain)
MIN_LOSS_BP = 100_000      # ... and at least this long
BREAK_WIN = 10_000         # evidence window at a loss's breakpoint
BREAK_WIN_MOSAIC = 20_000  # mosaic edges are refined to ~10 kb
TELO_WIN = 30_000          # telomeric clips this close: a depth-only edge is ~20 kb off
TELO_MIN = 3               # telomeric clips for a healed telomere
REPEAT_WIN = 2_000         # a junction to a repeat this close to the breakpoint is its own
TELO_REPEATS = 3           # TTAGGG / CCCTAA repeats in a clip to call it telomeric
INV_MIN_BP = 1_000         # shorter same-strand pairs are fold-back artefacts
INVDUP_SPACER_BP = 1_000_000  # an inv-dup's fold-back is this close to its loss
MAIN_CONTIGS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])


@dataclass
class Bnd:
    sample: str
    contig_a: str
    pos_a: int
    side_a: str           # "R": sequence left of pos_a is joined (clip on the right), "L", "?"
    contig_b: str
    pos_b: int
    side_b: str
    pe: int
    sr_a: int
    sr_b: int
    mapq: float
    flag: str = ""        # "", "unplaced", "subtelomeric", "promiscuous", "recurrent"
    kind: str = "TRA"     # "TRA" between contigs; same-strand pairs within one:
                          # "INV" (pairs on both sides of the breakpoints), "FOLD" (one side)


@dataclass
class TerminalLoss:
    sample: str
    contig: str
    start: int
    end: int
    source: str           # "call" / "mosaic"
    mechanism: str        # "translocation" / "inv-dup" / "telomere" / "repeat" / "unknown"
    partner: str          # chrB:pos of the junction, "" if none
    partner_gain: str     # the sample's gain at the partner's end (inv-dup: the
                          # inverted duplication), "" if none
    telo_clips: int
    junction_pe: int


_SIDE = {SR_RIGHT_CLIP: "R", SR_LEFT_CLIP: "L"}


class _Clips:
    """A sample's split-read pileups, made once per contig: a genome has
    ~13,000 pair clusters to place, those in repeats megabases wide, and
    piling up each one's clips anew was 90 % of the scan."""

    def __init__(self, store, config: BreakpointConfig):
        self.store, self.config, self._piles = store, config, {}

    def pileup(self, contig: str, lo: int, hi: int):
        """The strongest split-read pileup in [lo, hi] (the first of equals), or None."""
        if contig not in self._piles:
            bps = []
            if self.store.has_contig(contig):
                sr = self.store.sr(contig, 0, 2**32 - 1)
                if len(sr["pos"]):
                    bps = cluster_sr(sr, self.config)
            self._piles[contig] = (np.array([b.pos for b in bps], dtype=np.int64),
                                   np.array([b.n for b in bps], dtype=np.int64), bps)
        pos, n, bps = self._piles[contig]
        i, j = int(np.searchsorted(pos, max(lo, 0), "left")), int(np.searchsorted(pos, hi, "right"))
        if j <= i:
            return None
        return bps[i + int(np.argmax(n[i:j]))]


def scan_junctions(store, names: list[str], config: BreakpointConfig | None = None) -> list[Bnd]:
    """Inter-chromosomal and inversion junctions of one sample. `names`: the
    reference's contig names in index order (pairs store their mate's index)."""
    config = config or BreakpointConfig()
    cfg = BreakpointConfig(**{**asdict(config), "min_pe": BND_MIN_PE})
    insert = store.insert_model()
    mu = int(insert.mu) or 500
    clips = _Clips(store, config)
    found: dict[tuple, Bnd] = {}
    for contig in store.contigs:
        g = store.root[contig]
        if "pe" not in g:
            continue
        pe = g["pe"]
        orient = np.asarray(pe["orient"][:])
        found.update(_inversions(store, contig, pe, orient, insert, cfg, clips))
        sel = orient == PE_INTER
        if sel.sum() < BND_MIN_PE:
            continue
        cols = {k: np.asarray(pe[k][:])[sel] for k in ("pos1", "pos2", "mate_tid", "orient", "mapq")}
        tids, counts = np.unique(cols["mate_tid"], return_counts=True)
        for t in tids[counts >= BND_MIN_PE]:
            if not 0 <= t < len(names):
                continue
            mate = names[t]
            m = cols["mate_tid"] == t
            for cl in cluster_pe({k: v[m] for k, v in cols.items()}, insert, cfg, PE_INTER):
                if cl.mean_mapq < BND_MIN_MAPQ:
                    continue
                a = clips.pileup(contig, cl.pos1_lo - mu, cl.pos1_hi + mu)
                b = clips.pileup(mate, cl.pos2_lo - mu, cl.pos2_hi + mu)
                bnd = Bnd(
                    sample=store.sample, contig_a=contig,
                    pos_a=a.pos if a else (cl.pos1_lo + cl.pos1_hi) // 2,
                    side_a=_SIDE[a.side] if a else "?",
                    contig_b=mate,
                    pos_b=b.pos if b else (cl.pos2_lo + cl.pos2_hi) // 2,
                    side_b=_SIDE[b.side] if b else "?",
                    pe=cl.n, sr_a=a.n if a else 0, sr_b=b.n if b else 0,
                    mapq=round(cl.mean_mapq, 1),
                )
                key = (contig, bnd.pos_a, mate, bnd.pos_b)
                if key not in found or found[key].pe < bnd.pe:
                    found[key] = bnd   # two pair clusters on the same pileups
    return list(found.values())


def _inversions(store, contig, pe, orient, insert, cfg, clips) -> dict[tuple, Bnd]:
    """Same-strand pair clusters of one contig. Pairs are stored without
    strands: a junction's ++ pairs lie left of its first breakpoint, its --
    pairs right of it, so a reciprocal inversion has pairs on both sides and a
    fold-back (one junction of an inverted duplication) on one."""
    pos1 = np.asarray(pe["pos1"][:])
    pos2 = np.asarray(pe["pos2"][:])
    sel = (orient == PE_INV) & (pos2 - pos1 >= INV_MIN_BP)
    if sel.sum() < BND_MIN_PE:
        return {}
    cols = {"pos1": pos1[sel], "pos2": pos2[sel], "orient": orient[sel],
            "mapq": np.asarray(pe["mapq"][:])[sel]}
    mu = int(insert.mu) or 500
    out = {}
    for cl in cluster_pe(cols, insert, cfg, PE_INV):
        if cl.mean_mapq < BND_MIN_MAPQ:
            continue
        a = clips.pileup(contig, cl.pos1_lo - mu, cl.pos1_hi + mu)
        b = clips.pileup(contig, cl.pos2_lo - mu, cl.pos2_hi + mu)
        start = a.pos if a else (cl.pos1_lo + cl.pos1_hi) // 2
        m = ((cols["pos1"] >= cl.pos1_lo) & (cols["pos1"] <= cl.pos1_hi)
             & (cols["pos2"] >= cl.pos2_lo) & (cols["pos2"] <= cl.pos2_hi))
        left = int((cols["pos1"][m] < start).sum())
        both = min(left, int(m.sum()) - left) >= 2
        bnd = Bnd(
            sample=store.sample, contig_a=contig, pos_a=start,
            side_a=_SIDE[a.side] if a else "?", contig_b=contig,
            pos_b=b.pos if b else (cl.pos2_lo + cl.pos2_hi) // 2,
            side_b=_SIDE[b.side] if b else "?",
            pe=cl.n, sr_a=a.n if a else 0, sr_b=b.n if b else 0,
            mapq=round(cl.mean_mapq, 1), kind="INV" if both else "FOLD",
        )
        key = (contig, bnd.pos_a, contig, bnd.pos_b)
        if key not in out or out[key].pe < bnd.pe:
            out[key] = bnd
    return out


def flag_junctions(by_sample: dict[str, list[Bnd]], lengths: dict[str, int]) -> None:
    """Set each junction's flag from the run's genomes (see the module doc)."""
    def loci(b):
        return (b.contig_a, b.pos_a // LOCUS_BP), (b.contig_b, b.pos_b // LOCUS_BP)

    partners: dict[tuple[str, int], set] = defaultdict(set)
    pairs: dict[tuple, set[str]] = defaultdict(set)
    for s, bnds in by_sample.items():
        for b in bnds:
            a, z = loci(b)
            partners[a].add(z)
            partners[z].add(a)
            pairs[(a, z)].add(s)
    limit = max(BND_RECURRENT_MIN, BND_RECURRENT_FRAC * len(by_sample))

    def subtel(c, p):
        n = lengths.get(c)
        return n is not None and (p < SUBTEL_BP or p > n - SUBTEL_BP)

    for bnds in by_sample.values():
        for b in bnds:
            a, z = loci(b)
            if b.contig_a not in MAIN_CONTIGS or b.contig_b not in MAIN_CONTIGS:
                b.flag = "unplaced"
            elif subtel(b.contig_a, b.pos_a) or subtel(b.contig_b, b.pos_b):
                b.flag = "subtelomeric"
            elif max(len(partners[a]), len(partners[z])) >= BND_PROMISCUOUS_PARTNERS:
                b.flag = "promiscuous"
            elif len(pairs[(a, z)]) >= limit:
                b.flag = "recurrent"
            else:
                b.flag = ""


def telomeric(seq: bytes) -> bool:
    return seq.count(b"TTAGGG") >= TELO_REPEATS or seq.count(b"CCCTAA") >= TELO_REPEATS


def telomeric_clips(store, contig: str, pos: int, win: int, side: int) -> int:
    """Clipped reads within `win` of `pos`, clipped on `side`, made of telomere repeats."""
    sr = store.sr(contig, max(pos - win, 0), pos + win, with_seq=True)
    return sum(1 for sd, seq in zip(sr["side"], sr["seq"], strict=True) if int(sd) == side and telomeric(seq))


def callable_bins(store, contig: str) -> np.ndarray:
    """Start positions of the store's unmasked depth bins on a contig."""
    _, _, mask = store.depth(contig)
    return np.nonzero(~np.asarray(mask, dtype=bool))[0] * store.bin_size


def _chain(c: str, s: int, e: int, pieces: list[tuple], qter: bool) -> tuple[int, int]:
    """Extend a terminal loss inward over the abutting deletion pieces."""
    grown = True
    while grown:
        grown = False
        for pc, ps, pe in pieces:
            if pc != c or pe - ps < CHAIN_MIN_BP:
                continue
            if qter and ps < s and abs(pe - s) <= CHAIN_GAP_BP:
                s, grown = ps, True
            elif not qter and pe > e and abs(ps - e) <= CHAIN_GAP_BP:
                e, grown = pe, True
    return s, e


def terminal_losses(calls: list[tuple], mosaic: list[dict], lengths: dict[str, int],
                    callable_at=None, bin_size: int = 1000, pieces: list[tuple] = ()) -> list[tuple]:
    """(contig, start, end, source) of the losses reaching a contig end.
    `calls`: (contig, start, end) of PASS deletions; `mosaic`: mosaic rows.
    `pieces`: (contig, start, end) of all deletions, any filter: a large loss
    is often called in abutting pieces (a terminal 22q13 loss from 43.6 Mb
    as 43.6-49.3 + 49.3-50.8), and its breakpoint is that of the inner one.
    `callable_at(contig)`: its callable bin starts; with it, a loss is
    terminal only when the callable sequence past it is TERMINAL_TAIL_BP at
    most (the common 2q37.3, 6p25.3 and 7q36.3 deletions end 50-130 kb short)."""
    rows = [(c, s, e, "call") for c, s, e in calls]
    rows += [(m["contig"], int(m["start"]), int(m["end"]), "mosaic")
             for m in mosaic if m["svtype"] == "DEL"]
    out = []
    for c, s, e, src in rows:
        n = lengths.get(c)
        if n is None or not (e >= n - TERMINAL_BP or s <= TERMINAL_BP):
            continue
        if src == "call":
            s, e = _chain(c, s, e, [*calls, *pieces], e >= n - TERMINAL_BP)
        if e - s < MIN_LOSS_BP:
            continue
        if callable_at is not None:
            pos = callable_at(c)
            tail = (pos >= e).sum() if e >= n - TERMINAL_BP else (pos + bin_size <= s).sum()
            if tail * bin_size > TERMINAL_TAIL_BP:
                continue
        if (c, s, e, src) not in out:   # the pieces of one loss
            out.append((c, s, e, src))
    return out


def loss_mechanism(store, sample: str, loss: tuple, bnds: list[Bnd],
                   gains: list[tuple], lengths: dict[str, int]) -> TerminalLoss:
    """Mechanism of one terminal loss (see the module doc). `gains`:
    (contig, start, end) of the sample's gains (PASS duplications and mosaic
    gain segments)."""
    contig, start, end, source = loss
    qter = end >= lengths[contig] - TERMINAL_BP
    bp = start if qter else end
    # the retained side: left of a q-ter loss, so its reads clip on the right
    side = SR_RIGHT_CLIP if qter else SR_LEFT_CLIP
    win = BREAK_WIN_MOSAIC if source == "mosaic" else BREAK_WIN
    telo = telomeric_clips(store, contig, bp, TELO_WIN, side)
    near = []   # (junction, partner contig, partner pos, distance) by the end closer to bp
    for b in bnds:
        ends = [(abs(p - bp), oc, op) for c, p, oc, op in
                ((b.contig_a, b.pos_a, b.contig_b, b.pos_b), (b.contig_b, b.pos_b, b.contig_a, b.pos_a))
                if c == contig]
        if ends and min(ends)[0] <= win:
            dist, oc, op = min(ends)
            near.append((b, oc, op, dist))
    # a fold-back on the retained side: the inverted duplication of an inv-dup-del
    fold = sorted((x for x in near if x[0].kind != "TRA" and x[0].flag == ""
                   and 0 < (bp - x[2] if qter else x[2] - bp) <= INVDUP_SPACER_BP),
                  key=lambda x: -x[0].pe)
    trans = sorted((x for x in near if x[0].kind == "TRA" and x[0].flag == ""),
                   key=lambda x: -x[0].pe)
    subt = [x for x in near if x[0].kind == "TRA" and x[0].flag == "subtelomeric"]
    # a junction to a repeat at the breakpoint itself: a translocation onto
    # sequence the reference does not place (rDNA, acrocentric short arms)
    rep = sorted((x for x in near if x[0].kind == "TRA" and x[3] <= REPEAT_WIN
                  and x[0].flag in ("unplaced", "promiscuous")), key=lambda x: -x[0].pe)
    partner = gain = ""
    if trans:
        b, pc, pp, _ = trans[0]
        mechanism, partner, pe = "translocation", f"{pc}:{pp}", b.pe
        for gc, gs, ge in gains:
            if gc == pc and min(abs(gs - pp), abs(ge - pp)) <= 100_000:
                gain = f"{gc}:{gs}-{ge}"
    elif fold:
        b, pc, pp, _ = fold[0]
        mechanism, partner, pe = "inv-dup", f"{pc}:{pp}", b.pe
        for gc, gs, ge in gains:
            if gc == pc and abs((ge if qter else gs) - pp) <= 100_000:
                gain = f"{gc}:{gs}-{ge}"
    elif telo >= TELO_MIN or subt:
        mechanism, pe = "telomere", max((x[0].pe for x in subt), default=0)
    elif rep:
        b, pc, pp, _ = rep[0]
        mechanism, partner, pe = "repeat", f"{pc}:{pp}", b.pe
    else:
        mechanism, pe = "unknown", 0
    return TerminalLoss(sample=sample, contig=contig, start=start, end=end, source=source,
                        mechanism=mechanism, partner=partner, partner_gain=gain,
                        telo_clips=telo, junction_pe=pe)


BND_FIELDS = list(Bnd.__dataclass_fields__)
TERMINAL_FIELDS = list(TerminalLoss.__dataclass_fields__)


def write_tsv(rows: list, fields: list[str], out: str | Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def read_tsv(path: str | Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


MECHANISM_FR = {"translocation": "translocation", "inv-dup": "duplication inversée (inv-dup)",
                "telomere": "télomère néoformé",
                "repeat": "jonction vers une séquence répétée", "unknown": "mécanisme indéterminé"}


def read_terminal(path: str | Path) -> list[dict]:
    """A <S>.terminal.tsv, numbers as numbers."""
    rows = read_tsv(path)
    for r in rows:
        for k in ("start", "end", "telo_clips", "junction_pe"):
            r[k] = int(r[k])
    return rows


def _loci(contig_a, pos_a, contig_b, pos_b) -> frozenset:
    return frozenset({(contig_a, int(pos_a) // LOCUS_BP), (contig_b, int(pos_b) // LOCUS_BP)})


def loss_junction(row: dict, lengths: dict[str, int] | None = None) -> frozenset | None:
    """The two loci a terminal loss's junction joins (translocation, inv-dup), or None."""
    if not row["partner"]:
        return None
    pc, pp = row["partner"].rsplit(":", 1)
    n = (lengths or {}).get(row["contig"])
    qter = row["start"] > 0 if n is None else row["end"] >= n - TERMINAL_BP
    return _loci(row["contig"], row["start"] if qter else row["end"], pc, pp)


def junction_carriers(row: dict, others: dict[str, list[dict]],
                      lengths: dict[str, int] | None = None) -> list[str]:
    """The samples of `others` ({sample: <S>.bnd.tsv rows}) carrying the loss's
    junction: a parent with it carries the balanced form of the translocation
    (or the same inverted duplication)."""
    j = loss_junction(row, lengths)
    if j is None:
        return []
    return [s for s, bnds in others.items()
            if any(_loci(b["contig_a"], b["pos_a"], b["contig_b"], b["pos_b"]) == j for b in bnds)]


def loss_carriers(row: dict, others: dict[str, list[dict]]) -> list[str]:
    """The samples of `others` ({sample: <S>.terminal.tsv rows}) with the same
    terminal loss (same contig, reciprocal overlap >= 0.5): a parent with it
    transmitted it. The only two losses with identical bounds in two genomes
    of PMS and dataset_paper are a mother's and her child's."""
    def same(u):
        ov = min(u["end"], row["end"]) - max(u["start"], row["start"])
        return (u["contig"] == row["contig"]
                and ov >= 0.5 * max(u["end"] - u["start"], row["end"] - row["start"]))
    return [s for s, rows in others.items() if any(same(u) for u in rows)]


def terminal_label(row: dict, carriers: list[str] = (), inherited: list[str] = ()) -> str:
    """One line for a report: mechanism, partner, partner gain, carrier parent,
    parent with the same loss."""
    out = MECHANISM_FR.get(row["mechanism"], row["mechanism"])
    if row["partner"]:
        pc, pp = row["partner"].rsplit(":", 1)
        out += f" {pc}:{int(pp):,}"
    if row["partner_gain"]:
        gc, span = row["partner_gain"].rsplit(":", 1)
        gs, ge = (int(x) for x in span.split("-"))
        out += f", gain {gc}:{gs / 1e6:.1f}-{ge / 1e6:.1f} Mb"
    if carriers:
        out += f", jonction chez {', '.join(carriers)}"
    if inherited:
        out += f", héritée de {', '.join(inherited)}"
    return out


def _pass_calls(vcf: Path) -> tuple[list[tuple], list[tuple], list[tuple]]:
    """(losses, gains, pieces) as (contig, start, end): a filtered VCF's PASS
    deletions and duplications, and its deletions of any filter."""
    import pysam

    losses, gains, pieces = [], [], []
    if not vcf.exists():
        return losses, gains, pieces
    with pysam.VariantFile(str(vcf)) as v:
        for r in v:
            t = r.info.get("SVTYPE")
            if t == "DEL":
                pieces.append((r.chrom, r.start, r.stop))
            if list(r.filter.keys()) != ["PASS"]:
                continue
            (losses if t == "DEL" else gains if t == "DUP" else []).append((r.chrom, r.start, r.stop))
    return losses, gains, pieces


def read_bnds(path: str | Path) -> list[Bnd]:
    types = {k: f.type for k, f in Bnd.__dataclass_fields__.items()}
    conv = {"int": int, "float": float, "str": str}
    return [Bnd(**{k: conv[types[k]](v) for k, v in r.items()}) for r in read_tsv(path)]


def _scan_one(args):
    """A sample's junctions, unflagged; cached in <S>.junctions.tsv (the scan
    is ~1 min per genome, flagging needs the whole run)."""
    from .extract.store import EvidenceStore

    d, names = args
    store_path = d / f"{d.name}.evidence.zarr"
    cache = d / f"{d.name}.junctions.tsv"
    if cache.exists() and cache.stat().st_mtime >= store_path.stat().st_mtime:
        with open(cache) as fh:
            if fh.readline().rstrip("\n").split("\t") == BND_FIELDS:
                return read_bnds(cache)
    bnds = scan_junctions(EvidenceStore(store_path), names)
    write_tsv(bnds, BND_FIELDS, cache)
    return bnds


def _annotate_one(args):
    from .extract.store import EvidenceStore
    from .mosaic import read_tsv as read_mosaic
    from .vcfout import write_bnd_vcf

    d, bnds, lengths, contigs, fasta = args
    s = d.name
    store = EvidenceStore(d / f"{s}.evidence.zarr")
    losses, gains, pieces = _pass_calls(d / f"{s}.filtered.vcf.gz")
    mpath = d / f"{s}.mosaic.tsv"
    mosaic = read_mosaic(mpath) if mpath.exists() else []
    gains = gains + [(m["contig"], m["start"], m["end"]) for m in mosaic if m["svtype"] == "DUP"]
    rows = [loss_mechanism(store, s, loss, bnds, gains, lengths)
            for loss in terminal_losses(losses, mosaic, lengths, partial(callable_bins, store),
                                             store.bin_size, pieces)]
    # ~13,000 clusters per genome, nearly all in repeats: the files keep the rest
    bnds = [b for b in bnds if b.flag in ("", "subtelomeric")]
    write_tsv(bnds, BND_FIELDS, d / f"{s}.bnd.tsv")
    write_bnd_vcf(bnds, contigs, s, d / f"{s}.bnd.vcf.gz", fasta)
    write_tsv(rows, TERMINAL_FIELDS, d / f"{s}.terminal.tsv")
    return s, len(bnds), rows


def run_bnd(run_dir: str | Path, fasta: str | Path, jobs: int = 1) -> tuple[int, int, list]:
    """Every samples/<S>/ of a run: <S>.bnd.tsv (junctions, flagged over the
    run), <S>.bnd.vcf.gz (those with split reads at both ends, as breakend
    pairs) and <S>.terminal.tsv (terminal losses and their mechanism)."""
    from multiprocessing import get_context

    from .reference import read_fai

    contigs = read_fai(fasta)
    names = [c.name for c in contigs]
    lengths = {c.name: c.length for c in contigs}
    dirs = [d for d in sorted((Path(run_dir) / "samples").iterdir())
            if (d / f"{d.name}.evidence.zarr").exists()]
    with get_context("fork").Pool(max(jobs, 1)) as pool:
        found = pool.map(_scan_one, [(d, names) for d in dirs], chunksize=4)
        by_sample = {d.name: b for d, b in zip(dirs, found, strict=True)}
        flag_junctions(by_sample, lengths)
        res = pool.map(_annotate_one, [(d, by_sample[d.name], lengths, contigs, fasta)
                                       for d in dirs], chunksize=4)
    terminal = [r for _, _, rows in res for r in rows]
    return len(dirs), sum(n for _, n, _ in res), terminal
