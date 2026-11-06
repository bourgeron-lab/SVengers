"""Mosaic chromosomal-scale CNVs: fractional copy number, confirmed by BAF.

The CN-HMM has integer states only. A copy ratio between two of them — a loss
or gain present in a fraction of the cells — makes it alternate between the
neighbouring states: in a dataset_paper parent, a 50 Mb 10q terminal loss at
CN ~1.33 (two thirds of the cells) came back as 57 CN=1 fragments, and the
gains of the same clone (CN ~2.67 on a whole chromosome and a chromosome arm)
were not called at all.

This pass works at chromosomal scale, from the evidence store alone:

1. the diploid-normalised depth ratio is summarised in 1 Mb windows;
2. each autosome is segmented by binary segmentation on the window levels;
3. a segment away from CN 2 gets a cell fraction, assuming a single-copy
   change: loss CN = 2 - f, gain CN = 2 + f;
4. it must be homogeneous (HOMOG_*): a blend of integer levels is not a clone;
5. it is kept only if the B-allele frequency of its heterozygous sites departs
   from 0.5 as that fraction predicts — loss |BAF - 0.5| = 0.5 - (1-f)/(2-f),
   gain (1+f)/(2+f) - 0.5. GC or wave artefacts and pericentromeric repeats
   shift depth without any allelic imbalance, and fail here;
6. an edge inside the contig is refined to the store's bin (`refine_edge`).

The BAF statistic is the median |BAF - 0.5| of a window's heterozygous sites.
Two things move it besides a clone (`baf_rise`):
- per-site sampling noise (sd ~0.08 at 30x) absorbs most of a small shift:
  at f = 0.2 the theoretical shift is 0.056 but the median moves by ~0.013;
- depth itself: sampling noise grows as 1/sqrt(depth), so a region read at
  92 % depth (a GC wave) has a wider BAF spread and a median ~0.0025 higher
  with no allelic imbalance at all.
Both expectations are computed on the sample's own heterozygous sites — null:
noise scaled by the segment's depth; mosaic: the same plus ±(theoretical
shift) — and the observed excess must cover BAF_SUPPORT of the way from the
null to the mosaic expectation, BAF_Z window-noise SEs above the null. Only a
shallower segment widens the noise: not all of it is read sampling, and
crediting a deeper one with a narrower spread let GC-rich gain waves (BAF
excess exactly 0) through. Like depth, the BAF shift of a clone is the same
in every window (BAF_HOMOG_*): a segment whose median rides on one or two
LOH-like windows (excess 0.05-0.3) is not one.

On a wavy genome a segment needs a fraction of at least WAVY_MIN_FRAC.
Depth waves of ±15 % read as 15-25 % clones, and on GC-rich chromosomes such
a genome's heterozygous sites spread wider: once the BAF was read without a
QUAL threshold (`vedette baf`, twice the sites, a tighter window noise), 9
segments at f = 0.15-0.26 passed the BAF test on genomes with a wave score of
0.134-0.210 (chr16, chr17, chr19, chr22, 7q11.23), the carrier's own excess
in all but chr16, where a 200-genome panel is as high. The mosaics known
before sit at f >= 0.38 on wavy genomes, f = 0.16-0.21 on genomes below 0.1.

A fixed floor of 0.03 missed a PMS proband's 20 % 22q13 loss (41.7 Mb to the
end, observed 0.012, null 0.003, mosaic 0.016, seen in IGV); without the
depth null, dataset_paper gained ~70 wave segments at CN ~1.83 (chr4, chr13,
chr19 ...), whose excess (0.003-0.006) is the null's.

Autosomes only: the sex chromosomes need the sample's ploidy, which older
stores do not record.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .qc import WAVE_WARN

WINDOW_BP = 1_000_000
# 3 Mb. At 5 Mb a PMS proband's terminal 22q13 loss at CN 1.62 (47.7 Mb to
# the end of chr22, BAF exactly as f = 0.38 predicts, confirmed in IGV) was
# out of reach; at 3 windows the PMS scan finds three terminal 22q13 mosaics
# and nothing else.
MIN_WINDOWS = 3
MIN_VALID_FRAC = 0.5       # a window needs half its bins unmasked
MIN_HET_PER_WINDOW = 20    # heterozygous sites for a window's BAF summary
SPLIT_T = 6.0              # binary-segmentation split threshold (in noise SDs)
MIN_CN_SHIFT = 0.15        # |CN - 2| below this is "no change"
MOSAIC_FRAC = (0.15, 0.85)
BAF_SUPPORT = 0.5          # share of the way from the depth null to the mosaic expectation
BAF_Z = 5.0                # ... and this many window-noise SEs above the depth null
BAF_HOMOG_MAX = 3.0        # a window's excess beyond 3x the expected is not this clone
BAF_HOMOG_FRAC = 0.5       # share of the windows between the floor and that ceiling
# A mosaic is a homogeneous intermediate level. A segment that blends windows
# at different integer levels has the same median and fails here: on PMS, a
# 3 Mb constitutional 16p13 duplication next to a deleted window (2.9, 2.5,
# 3.0, 1.15) and a 22q13 deletion with an intermediate edge (1.66, 1.65,
# 1.03, 1.02) both read as "mosaic" without it.
HOMOG_TOL = 0.2
HOMOG_FRAC = 0.75
WAVY_MIN_FRAC = 0.3        # least cell fraction on a wavy genome (see above)

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


@dataclass
class MosaicCall:
    sample: str
    contig: str
    start: int            # 0-based; an inner edge is refined to the store bin
    end: int
    svtype: str           # DEL (loss) / DUP (gain)
    cn: float
    fraction: float       # estimated cell fraction
    baf_dev: float        # observed median |BAF - 0.5| minus the genome baseline
    baf_expected: float
    n_windows: int
    mosaic: bool          # fraction within MOSAIC_FRAC

    @property
    def size(self) -> int:
        return self.end - self.start


def expected_baf_excess(svtype: str, f: float) -> float:
    if svtype == "DEL":
        return 0.5 - (1 - f) / (2 - f)
    return (1 + f) / (2 + f) - 0.5


def baf_rise(het_dev: np.ndarray, depth_ratio: float, shift: float = 0.0) -> float:
    """Expected rise of the median |BAF - 0.5| in a segment read at
    `depth_ratio` of the genome's depth whose heterozygous sites are shifted
    by ±`shift`, from the sample's own deviations `het_dev` (BAF - 0.5)."""
    if len(het_dev) == 0:
        return shift
    noise = het_dev / np.sqrt(min(max(depth_ratio, 0.05), 1.0))
    shifted = np.concatenate([np.abs(noise + shift), np.abs(noise - shift)])
    return float(np.median(shifted) - np.median(np.abs(het_dev)))


def contig_wave(meta: dict, contig: str, min_contigs: int = 10) -> float | None:
    """The genome's wave score without ``contig`` (a large CNV raises its own
    contig's score), the sample score on stores without per-contig values."""
    per = meta.get("wave_by_contig") or {}
    vals = [v for k, v in per.items() if k != contig and v is not None]
    if len(vals) >= min_contigs:
        return float(np.mean(vals))
    w = meta.get("wave_score")
    return None if w is None else float(w)


def refine_edge(pos: np.ndarray, y: np.ndarray, edge: int, inner: float,
                bin_size: int, left: bool, min_step: float = 0.0) -> int:
    """Breakpoint of a segment edge found at window resolution, to the bin.

    The true edge lies within a window of `edge`: the segmenter cut where the
    window levels change, and the window that straddles the breakpoint goes to
    either side. Over the bins (`pos`, `y`: valid bins of the contig, depth
    not normalised) of the two windows around it, the step from the outer
    level (the median of the next window out) to the segment's own level
    `inner` is placed where the cumulative evidence for the outer level peaks;
    each bin counts within ±(step height), so a spike cannot move it. A 15 %
    clone changes 1 kb depth by ~0.08 against a bin noise of ~0.12: a PMS
    proband's 22q13 loss, scanned at 42.0 Mb, comes back at 41,993,000, 6 kb
    from its telomere-healing junction (one split read, 41,999,335).

    `left`: the segment lies right of `edge` (its start). The edge is left as
    it is when the next window out has no valid bins (a centromere, a gap),
    when that window is within `min_step` of the segment's level, or when the
    step lands on a limit of the search: the breakpoint is then out of reach
    (a PMS 22q13 duplication, 41.5-49.25 Mb, scanned at 44-48 Mb next to two
    blended segments the homogeneity test dropped)."""
    w = WINDOW_BP
    ref = (pos >= edge - 2 * w) & (pos < edge - w) if left else \
        (pos >= edge + w) & (pos < edge + 2 * w)
    near = (pos >= edge - w) & (pos < edge + w)
    if ref.sum() < 0.1 * w / bin_size or not near.any():
        return edge
    outer = float(np.median(y[ref]))
    step = abs(outer - inner)
    if step == 0 or step < min_step:
        return edge
    x, v = pos[near], y[near]
    # positive where a bin looks like the outer level
    c = np.clip(np.sign(outer - inner) * (v - (outer + inner) / 2), -step, step)
    if left:
        score = np.concatenate([[0.0], np.cumsum(c)])    # outer bins = [0, k)
    else:
        score = np.concatenate([np.cumsum(c[::-1])[::-1], [0.0]])  # outer bins = [k, n)
    k = int(np.argmax(score))
    if k == 0 or k == len(x):
        return edge
    return int(x[k]) if left else int(x[k - 1] + bin_size)


def homogeneous(win_cn: np.ndarray, cn: float) -> bool:
    """Whether a segment's windows sit at one level (HOMOG_*), as a clone does."""
    return bool(np.mean(np.abs(np.asarray(win_cn) - cn) <= HOMOG_TOL) >= HOMOG_FRAC)


def _segments(levels: np.ndarray, sd: float, min_len: int, t_split: float) -> list[tuple[int, int]]:
    """Binary segmentation of window levels into [i, j) runs."""
    out: list[tuple[int, int]] = []
    stack = [(0, len(levels))]
    while stack:
        i, j = stack.pop()
        x = levels[i:j]
        n = len(x)
        best_k, best_t = None, 0.0
        if n >= 2 * min_len:
            cs = np.cumsum(x)
            tot = cs[-1]
            for k in range(min_len, n - min_len + 1):
                ml = cs[k - 1] / k
                mr = (tot - cs[k - 1]) / (n - k)
                t = abs(ml - mr) / (sd * np.sqrt(1.0 / k + 1.0 / (n - k)))
                if t > best_t:
                    best_k, best_t = k, t
        if best_k is not None and best_t >= t_split:
            stack.append((i + best_k, j))
            stack.append((i, i + best_k))
        else:
            out.append((i, j))
    out.sort()
    return out


def _windows(store, contig: str):
    """(starts, ratio level, BAF dev, length, het BAF - 0.5) per 1 Mb window
    of one contig, ratio not yet normalised, and (pos, depth) of its valid
    bins for `refine_edge`."""
    _, corr, mask = store.depth(contig)
    corr = np.asarray(corr, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    nb = max(1, WINDOW_BP // store.bin_size)
    # keep the last, partial window: a terminal event lives in it
    n_win = -(-len(corr) // nb)
    pad = n_win * nb - len(corr)
    corr = np.concatenate([corr, np.full(pad, np.nan)])
    mask = np.concatenate([mask, np.ones(pad, dtype=bool)])
    blk = corr.reshape(n_win, nb).copy()
    bad = mask.reshape(n_win, nb) | ~np.isfinite(blk)
    blk[bad] = np.nan
    valid = (~bad).mean(axis=1) >= MIN_VALID_FRAC
    bins = (np.nonzero(~bad.ravel())[0] * store.bin_size, blk.ravel()[~bad.ravel()])
    level = np.full(n_win, np.nan)
    if valid.any():
        level[valid] = np.nanmedian(blk[valid], axis=1)
    dev = np.full(n_win, np.nan)
    het_dev = np.empty(0)
    try:
        baf = store.baf(contig, 0, n_win * WINDOW_BP)
        pos = np.asarray(baf["pos"], dtype=np.int64)
        vaf = np.asarray(baf["vaf"], dtype=np.float64)
        het = (vaf > 0.1) & (vaf < 0.9)
        pos, het_dev = pos[het], vaf[het] - 0.5
        a = np.abs(het_dev)
        order = np.argsort(pos, kind="stable")
        pos, a = pos[order], a[order]
        cuts = np.searchsorted(pos, np.arange(n_win + 1) * WINDOW_BP)
        for w in range(n_win):
            if cuts[w + 1] - cuts[w] >= MIN_HET_PER_WINDOW:
                dev[w] = np.median(a[cuts[w]:cuts[w + 1]])
    except KeyError:
        pass
    starts = np.arange(n_win) * WINDOW_BP
    return starts, level, dev, (len(corr) - pad) * store.bin_size, het_dev, bins


def scan_store(store) -> list[MosaicCall]:
    """Chromosomal-scale CNV segments of one sample, with their cell fraction
    and BAF support; only BAF-supported segments are returned."""
    contigs = [c for c in AUTOSOMES if store.has_contig(c)]
    data = {c: _windows(store, c) for c in contigs}
    length = {c: d[3] for c, d in data.items()}
    all_levels = np.concatenate([d[1] for d in data.values()]) if data else np.array([])
    base = np.nanmedian(all_levels) if np.isfinite(all_levels).any() else np.nan
    if not np.isfinite(base) or base <= 0:
        return []
    all_dev = np.concatenate([d[2] for d in data.values()])
    fin_dev = all_dev[np.isfinite(all_dev)]
    baf_base = np.median(fin_dev) if len(fin_dev) else np.nan
    baf_sd = max(1.4826 * np.median(np.abs(fin_dev - baf_base)), 1e-3) if len(fin_dev) else np.nan
    het_dev = np.concatenate([d[4] for d in data.values()])
    # window-level noise: MAD of successive differences, robust to the steps
    diffs = np.concatenate([np.diff(d[1][np.isfinite(d[1])] / base) for d in data.values()])
    sd = 1.4826 * np.median(np.abs(diffs - np.median(diffs))) / np.sqrt(2) if len(diffs) else 0.0
    sd = max(sd, 0.01)
    calls: list[MosaicCall] = []
    for c, (starts, level, dev, _len, _het, (bpos, bdep)) in data.items():
        ok = np.isfinite(level)
        if ok.sum() < MIN_WINDOWS:
            continue
        idx = np.nonzero(ok)[0]
        lv = level[ok] / base
        for i, j in _segments(lv, sd, MIN_WINDOWS, SPLIT_T):
            cn = 2.0 * float(np.median(lv[i:j]))
            if abs(cn - 2.0) < MIN_CN_SHIFT:
                continue
            if not homogeneous(2.0 * lv[i:j], cn):
                continue                   # blend of levels, not one clone
            svtype = "DEL" if cn < 2 else "DUP"
            f = float(np.clip(abs(cn - 2.0), 0.0, 1.0))
            w = contig_wave(store.meta, c)
            if f < WAVY_MIN_FRAC and w is not None and w >= WAVE_WARN:
                continue                   # a wave, not a clone
            wins = idx[i:j]
            d = dev[wins]
            if not np.isfinite(baf_base) or np.isfinite(d).sum() < 0.5 * len(wins):
                continue                   # no allelic evidence (repeats, gaps)
            excess = float(np.nanmedian(d) - baf_base)
            null = baf_rise(het_dev, cn / 2.0)
            exp = baf_rise(het_dev, cn / 2.0, expected_baf_excess(svtype, f))
            lo = null + BAF_SUPPORT * (exp - null)
            hi = null + BAF_HOMOG_MAX * (exp - null)
            wd = d[np.isfinite(d)] - baf_base
            z = (excess - null) / (baf_sd / np.sqrt(len(wd)))
            if excess < lo or z < BAF_Z or np.mean((wd >= lo) & (wd <= hi)) < BAF_HOMOG_FRAC:
                continue
            start = int(starts[wins[0]])
            end = int(min(starts[wins[-1]] + WINDOW_BP, length[c]))
            inner = cn / 2.0 * base
            # the window out must differ from the segment: a duplication that
            # runs on past the window (two PMS genomes) has no step there.
            # A quarter of the shift: CN 4 next to CN 3 is half of it.
            min_step = 0.25 * abs(inner - base)
            if start > 0:
                start = refine_edge(bpos, bdep, start, inner, store.bin_size, True, min_step)
            if end < length[c]:
                end = min(refine_edge(bpos, bdep, end, inner, store.bin_size, False, min_step),
                          length[c])
            calls.append(MosaicCall(
                sample=store.sample, contig=c,
                start=start, end=end,
                svtype=svtype, cn=round(cn, 3), fraction=round(f, 3),
                baf_dev=round(excess, 4), baf_expected=round(exp, 4),
                n_windows=int(len(wins)),
                mosaic=MOSAIC_FRAC[0] <= f <= MOSAIC_FRAC[1],
            ))
    return calls


FIELDS = ["sample", "contig", "start", "end", "svtype", "cn", "fraction",
          "baf_dev", "baf_expected", "n_windows", "mosaic"]


def write_tsv(calls: list[MosaicCall], out: str | Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for call in calls:
            row = asdict(call)
            row["mosaic"] = "1" if call.mosaic else "0"
            w.writerow(row)


def read_tsv(path: str | Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["start"], r["end"], r["n_windows"] = int(r["start"]), int(r["end"]), int(r["n_windows"])
        for k in ("cn", "fraction", "baf_dev", "baf_expected"):
            r[k] = float(r[k])
        r["mosaic"] = r["mosaic"] == "1"
    return rows
