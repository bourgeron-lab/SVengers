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
   shift depth without any allelic imbalance, and fail here.

Autosomes only: the sex chromosomes need the sample's ploidy, which older
stores do not record.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

WINDOW_BP = 1_000_000
# 3 Mb. At 5 Mb a PMS proband's terminal 22q13 loss at CN 1.62 (47.7 Mb to
# the end of chr22, BAF exactly as f = 0.38 predicts, confirmed in IGV) was
# out of reach; at 3 windows the PMS scan finds three terminal 22q13 mosaics
# and nothing else.
MIN_WINDOWS = 3
MIN_VALID_FRAC = 0.5       # a window needs half its bins unmasked
MIN_HET_PER_WINDOW = 20    # heterozygous sites for a window's BAF summary
SPLIT_T = 6.0              # binary-segmentation split threshold (in noise SDs)
MIN_CN_SHIFT = 0.2         # |CN - 2| below this is "no change"
MOSAIC_FRAC = (0.15, 0.85)
BAF_SUPPORT = 0.4          # observed BAF excess must reach this share of the expected
BAF_MIN_EXCESS = 0.03
# A mosaic is a homogeneous intermediate level. A segment that blends windows
# at different integer levels has the same median and fails here: on PMS, a
# 3 Mb constitutional 16p13 duplication next to a deleted window (2.9, 2.5,
# 3.0, 1.15) and a 22q13 deletion with an intermediate edge (1.66, 1.65,
# 1.03, 1.02) both read as "mosaic" without it.
HOMOG_TOL = 0.2
HOMOG_FRAC = 0.75

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


@dataclass
class MosaicCall:
    sample: str
    contig: str
    start: int            # 0-based, window-aligned
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
    """(starts, ratio level, BAF dev, valid) per 1 Mb window of one contig,
    ratio not yet normalised."""
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
    level = np.full(n_win, np.nan)
    if valid.any():
        level[valid] = np.nanmedian(blk[valid], axis=1)
    dev = np.full(n_win, np.nan)
    try:
        baf = store.baf(contig, 0, n_win * WINDOW_BP)
        pos = np.asarray(baf["pos"], dtype=np.int64)
        vaf = np.asarray(baf["vaf"], dtype=np.float64)
        het = (vaf > 0.1) & (vaf < 0.9)
        pos, a = pos[het], np.abs(vaf[het] - 0.5)
        order = np.argsort(pos, kind="stable")
        pos, a = pos[order], a[order]
        cuts = np.searchsorted(pos, np.arange(n_win + 1) * WINDOW_BP)
        for w in range(n_win):
            if cuts[w + 1] - cuts[w] >= MIN_HET_PER_WINDOW:
                dev[w] = np.median(a[cuts[w]:cuts[w + 1]])
    except KeyError:
        pass
    starts = np.arange(n_win) * WINDOW_BP
    return starts, level, dev, (len(corr) - pad) * store.bin_size


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
    baf_base = np.nanmedian(all_dev) if np.isfinite(all_dev).any() else np.nan
    # window-level noise: MAD of successive differences, robust to the steps
    diffs = np.concatenate([np.diff(d[1][np.isfinite(d[1])] / base) for d in data.values()])
    sd = 1.4826 * np.median(np.abs(diffs - np.median(diffs))) / np.sqrt(2) if len(diffs) else 0.0
    sd = max(sd, 0.01)
    calls: list[MosaicCall] = []
    for c, (starts, level, dev, _len) in data.items():
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
            wins = idx[i:j]
            d = dev[wins]
            if not np.isfinite(baf_base) or np.isfinite(d).sum() < 0.5 * len(wins):
                continue                   # no allelic evidence (repeats, gaps)
            excess = float(np.nanmedian(d) - baf_base)
            exp = expected_baf_excess(svtype, f)
            if excess < max(BAF_MIN_EXCESS, BAF_SUPPORT * exp):
                continue
            calls.append(MosaicCall(
                sample=store.sample, contig=c,
                start=int(starts[wins[0]]),
                end=int(min(starts[wins[-1]] + WINDOW_BP, length[c])),
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
