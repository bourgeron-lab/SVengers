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

COLUMNS = [
    "run", "svtype", "size", "median_per_child", "p90_per_child", "total",
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

    def median(self, key: tuple[str, str]) -> float:
        return statistics.median(self.per_child[key]) if self.per_child[key] else 0.0

    def quantile(self, key: tuple[str, str], q: float) -> int:
        v = sorted(self.per_child[key])
        return v[int(q * (len(v) - 1))] if v else 0


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
    )


def _rise(base: float, cand: float, max_increase: float, min_abs: float) -> bool:
    return cand - base >= min_abs and cand > base * (1 + max_increase)


def gate(
    baseline: RunStats,
    candidate: RunStats,
    max_increase: float = 0.25,
    min_abs: float = 2.0,
) -> tuple[list[dict], list[str]]:
    """Table rows for ``candidate`` against ``baseline`` and the alarms raised."""
    rows, alarms = [], []
    for key in BINS:
        b, c = baseline.median(key), candidate.median(key)
        flag = _rise(b, c, max_increase, min_abs)
        if flag:
            alarms.append(f"de novo {key[0]} {key[1]}: median per child {b:g} -> {c:g}")
        rows.append({
            "run": candidate.name, "svtype": key[0], "size": key[1],
            "median_per_child": c, "p90_per_child": candidate.quantile(key, 0.9),
            "total": sum(candidate.per_child[key]), "baseline_median": b,
            "change": f"{c - b:+g}", "flag": "ALARM" if flag else "",
        })
    b, c = baseline.recurrent_sites, candidate.recurrent_sites
    flag = _rise(b, c, max_increase, 10 * min_abs)
    if flag:
        alarms.append(f"recurrent de novo sites {b} -> {c}")
    rows.append({
        "run": candidate.name, "svtype": "recurrent_sites", "size": "all",
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
    return out


def write_table(rows: list[dict], out: str | Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
