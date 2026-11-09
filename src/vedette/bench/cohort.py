"""Cohort-level gate for a model or pipeline change: de novo burden per child.

A filter model can look better on a truth-set holdout, on the per-genome PASS
median of a panel and on "child call with a parent call" trio shares, and
still pass a family-shared artifact class (GC waves, segdup depth shifts) in
every genome of a family: such calls look inherited at the call level. Once
the whole cohort is genotyped they come out as de novo in hundreds of
children and as recurrent de novo sites. This module compares cohort runs
(Nextflow layout: ``<run>/cohort/{cohort.vcf.gz,mendel_report.tsv,
mendel_summary.tsv}``) on that burden:

* de novo calls per child at PASS sites, child and both parents GQ >= min_gq,
  by SVTYPE and size class (median, p90 over every child of the summary);
* recurrent de novo sites (de novo in >= ``recurrent`` children) and the calls
  they hold;
* PASS sites by type, median Mendelian-violation and transmission rates.

``gate`` flags a candidate whose median de novo burden (overall, or in any
type x size class) or recurrent-site count rises by more than ``max_increase``
over the baseline AND by at least ``min_abs`` calls per child (sites for the
recurrence), so that a median moving 0 -> 1 is not an alarm.

The same rule is applied within strata of the children's wave score (read
from ``<run>/samples/<child>/<child>.evidence.zarr``): a model that stops
using the wave score can stay quiet on the cohort median while flooding the
few wavy children (a candidate took 27 children at wave >= 0.20 from 18 to
460 de novo DUPs >= 1 kb each), and the strata also show whether a rise is
confined to them. A stratum needs ``MIN_STRATUM`` children to be gated.
"""

from __future__ import annotations

import csv
import logging
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)

SIZE_CLASSES = ("<1kb", "1-10kb", ">=10kb")
SVTYPES = ("DEL", "DUP")
BINS = [("all", "all")] + [(t, s) for t in SVTYPES for s in ("all", *SIZE_CLASSES)]

# wave QC strata: below the calibrated warning threshold, wavy, strongly wavy
WAVE_STRATA = (("<0.128", 0.0, 0.128), ("0.128-0.20", 0.128, 0.20), (">=0.20", 0.20, float("inf")))
MIN_STRATUM = 10

COLUMNS = [
    "run", "wave", "n_children", "svtype", "size", "median_per_child", "p90_per_child", "total",
    "baseline_median", "change", "flag",
]


def size_class(length: int) -> str:
    return "<1kb" if length < 1_000 else "1-10kb" if length < 10_000 else ">=10kb"


def cohort_dir(path: str | Path) -> Path:
    """Accept a run directory (``<run>/cohort/``) or the cohort directory itself."""
    p = Path(path)
    return p / "cohort" if (p / "cohort" / "mendel_report.tsv").exists() else p


@dataclass
class RunStats:
    name: str
    n_children: int
    per_child: dict[tuple[str, str], list[int]]
    recurrent_sites: int
    recurrent_calls: int
    pass_sites: Counter = field(default_factory=Counter)
    mvr: float = float("nan")
    transmission: float = float("nan")
    children: list[str] = field(default_factory=list)
    wave: dict[str, float] = field(default_factory=dict)

    def values(self, key: tuple[str, str], only: set[str] | None = None) -> list[int]:
        if only is None:
            return self.per_child[key]
        return [v for c, v in zip(self.children, self.per_child[key], strict=True) if c in only]

    def median(self, key: tuple[str, str], only: set[str] | None = None) -> float:
        v = self.values(key, only)
        return statistics.median(v) if v else 0.0

    def quantile(self, key: tuple[str, str], q: float, only: set[str] | None = None) -> int:
        v = sorted(self.values(key, only))
        return v[int(q * (len(v) - 1))] if v else 0


def read_waves(run_dir: str | Path, samples: list[str]) -> dict[str, float]:
    """{sample: wave_score} from the run's per-sample evidence stores (absent ones skipped)."""
    import zarr

    out = {}
    for s in samples:
        p = Path(run_dir) / "samples" / s / f"{s}.evidence.zarr"
        if not p.exists():
            continue
        try:
            w = zarr.open_group(str(p), mode="r").attrs.get("wave_score")
        except Exception:  # noqa: BLE001 - a broken store only loses its stratum
            continue
        if w is not None:
            out[s] = float(w)
    return out


def wave_strata(waves: dict[str, float]) -> list[tuple[str, set[str]]]:
    return [(name, {c for c, w in waves.items() if lo <= w < hi}) for name, lo, hi in WAVE_STRATA]


def pass_sites(cohort_vcf: str | Path) -> tuple[set[str], Counter]:
    ids: set[str] = set()
    n: Counter = Counter()
    with pysam.VariantFile(str(cohort_vcf), drop_samples=True) as vcf:
        for rec in vcf:
            if list(rec.filter.keys()) == ["PASS"]:
                ids.add(rec.id)
                n[rec.info.get("SVTYPE", "?")] += 1
    return ids, n


def load_run(
    path: str | Path, name: str | None = None, min_gq: int = 20, recurrent: int = 5
) -> RunStats:
    d = cohort_dir(path)
    name = name or (Path(path).name if d != Path(path) else Path(path).parent.name)
    summary = list(csv.DictReader(open(d / "mendel_summary.tsv"), delimiter="\t"))
    children = [r["child"] for r in summary]
    ids, n_pass = pass_sites(d / "cohort.vcf.gz")

    counts: dict[str, Counter] = defaultdict(Counter)
    per_site: Counter = Counter()
    with open(d / "mendel_report.tsv") as fh:
        col = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        i_id, i_inh, i_child = col["site_id"], col["inheritance"], col["child"]
        i_s, i_e, i_t = col["start"], col["end"], col["svtype"]
        i_gq = [col["child_gq"], col["father_gq"], col["mother_gq"]]
        for line in fh:
            if "\tdenovo\t" not in line:  # fast path: the report is ~10 GB
                continue
            f = line.rstrip("\n").split("\t")
            if f[i_inh] != "denovo" or f[i_id] not in ids:
                continue
            if min(int(f[i] or 0) for i in i_gq) < min_gq:
                continue
            sc, t = size_class(int(f[i_e]) - int(f[i_s])), f[i_t]
            c = counts[f[i_child]]
            for key in (("all", "all"), (t, "all"), (t, sc)):
                c[key] += 1
            per_site[f[i_id]] += 1

    rec = [k for k in per_site.values() if k >= recurrent]
    num = lambda r, k: float(r[k]) if r.get(k) not in (None, "", "NA") else None  # noqa: E731
    mvr = [x for r in summary if (x := num(r, "violation_rate")) is not None]
    tr = [x for r in summary if (x := num(r, "transmission_rate")) is not None]
    return RunStats(
        name=name,
        n_children=len(children),
        per_child={k: [counts[c][k] for c in children] for k in BINS},
        recurrent_sites=len(rec),
        recurrent_calls=sum(rec),
        pass_sites=n_pass,
        mvr=statistics.median(mvr) if mvr else float("nan"),
        transmission=statistics.median(tr) if tr else float("nan"),
        children=children,
        wave=read_waves(d.parent, children),
    )


def _rise(base: float, cand: float, max_increase: float, min_abs: float) -> bool:
    return cand - base >= min_abs and cand > base * (1 + max_increase)


def gate(
    baseline: RunStats,
    candidate: RunStats,
    max_increase: float = 0.25,
    min_abs: float = 2.0,
) -> tuple[list[dict], list[str]]:
    """Table rows for ``candidate`` against ``baseline`` and the alarms raised:
    over every child, then within each wave stratum of >= MIN_STRATUM children
    present in both runs (the wave is a property of the genome: the
    baseline's stores fill in for a candidate without them)."""
    rows, alarms = [], []
    waves = {**baseline.wave, **candidate.wave}
    common = set(baseline.children) & set(candidate.children)
    groups: list[tuple[str, set[str] | None]] = [("all", None)]
    groups += [(n, m & common) for n, m in wave_strata(waves) if len(m & common) >= MIN_STRATUM]
    for wname, only in groups:
        n = len(common) if only is None else len(only)
        for key in BINS:
            b, c = baseline.median(key, only), candidate.median(key, only)
            flag = _rise(b, c, max_increase, min_abs)
            if flag:
                where = "" if only is None else f" [wave {wname}, {n} children]"
                alarms.append(f"de novo {key[0]} {key[1]}{where}: median per child {b:g} -> {c:g}")
            rows.append({
                "run": candidate.name, "wave": wname, "n_children": n,
                "svtype": key[0], "size": key[1],
                "median_per_child": c, "p90_per_child": candidate.quantile(key, 0.9, only),
                "total": sum(candidate.values(key, only)), "baseline_median": b,
                "change": f"{c - b:+g}", "flag": "ALARM" if flag else "",
            })
    b, c = baseline.recurrent_sites, candidate.recurrent_sites
    flag = _rise(b, c, max_increase, 10 * min_abs)
    if flag:
        alarms.append(f"recurrent de novo sites {b} -> {c}")
    rows.append({
        "run": candidate.name, "wave": "all", "n_children": len(common),
        "svtype": "recurrent_sites", "size": "all",
        "median_per_child": "", "p90_per_child": "", "total": c,
        "baseline_median": b, "change": f"{c - b:+d}", "flag": "ALARM" if flag else "",
    })
    return rows, alarms


def summary_lines(s: RunStats) -> list[str]:
    out = [
        f"== {s.name}: {s.n_children} children, PASS sites "
        + ", ".join(f"{t} {n}" for t, n in sorted(s.pass_sites.items())),
        f"   MVR median {s.mvr:.4f}, transmission median {s.transmission:.3f}, "
        f"recurrent de novo sites {s.recurrent_sites} ({s.recurrent_calls} calls)",
    ]
    for key in BINS:
        out.append(
            f"   de novo {key[0]:>3} {key[1]:>6}: median {s.median(key):g}, "
            f"p90 {s.quantile(key, 0.9)}, total {sum(s.per_child[key])}"
        )
    for name, only in wave_strata(s.wave):
        if not only:
            continue
        out.append(
            f"   wave {name:>10} ({len(only)} children): de novo median all "
            f"{s.median(('all', 'all'), only):g}, DEL {s.median(('DEL', 'all'), only):g}, "
            f"DUP {s.median(('DUP', 'all'), only):g} (1-10kb {s.median(('DUP', '1-10kb'), only):g}, "
            f">=10kb {s.median(('DUP', '>=10kb'), only):g})"
        )
    if len(s.wave) < s.n_children:
        out.append(f"   wave score missing for {s.n_children - len(s.wave)} children (no store)")
    return out


def write_table(rows: list[dict], out: str | Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
