"""Pre-gate: per-genome filter output of a candidate model against a baseline,
before the cohort run.

The cohort gate (``bench cohort``) needs a whole cohort merged, genotyped and
PED-annotated (~2.5 h for 2,670 genomes); the per-sample filter step alone
takes minutes, and the de novo floods it catches start there:

* calls a candidate newly PASSes nominate cohort sites, which are then
  genotyped in every child and come out de novo wherever the depth genotyper
  sees a carrier with reference parents;
* LowSQ-only calls at or above the merge's ``rescue_min_sq`` join clusters as
  sub-threshold members and can rescue a site on their own (RESCUE=RECURRENT),
  so a model that lifts the SQ of the calls it still rejects floods the cohort
  without passing them.

A candidate that took the DUP 1-10 kb PASS calls of 2,670 genomes from 812 to
5,482 and their rescuable LowSQ calls from 45k to 136k (split-read-only DUPs
labelled true by transmission in wavy families) raised the children's de novo
DUP 1-10 kb median from 1 to 7; the panel of 186 genomes had shown its DUPs
going DOWN by 11 % per genome. On the seven candidates the cohort gate had
judged on that cohort, this pre-gate gives the same verdict for all seven.

Both counts are compared as a mean per genome, by SVTYPE x size class, over
every genome common to both runs and within the wave strata of
``bench.cohort``. Inputs are run directories with the Nextflow layout
``<run>/samples/<S>/<S>.filtered.vcf.gz`` (a ``--steps filter`` run is enough).
"""

from __future__ import annotations

import csv
import statistics
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .cohort import MIN_STRATUM, SIZE_CLASSES, SVTYPES, _rise, read_waves, size_class, wave_strata

KEYS = [(t, s) for t in SVTYPES for s in ("all", *SIZE_CLASSES)]
KINDS = ("pass", "rescuable")

COLUMNS = [
    "run", "wave", "n_genomes", "kind", "svtype", "size", "mean_per_genome", "median_per_genome",
    "total", "baseline_mean", "change", "flag",
]


def count_sample(vcf_path: str | Path, rescue_min_sq: float = 1.0) -> dict[str, Counter]:
    """{"pass": Counter, "rescuable": Counter} keyed by (svtype, size class) and
    (svtype, "all"), over the sample's carrier DEL/DUP records: PASS, and
    LowSQ-only with SQ >= rescue_min_sq (what the merge loads as sub-threshold
    members)."""
    out = {k: Counter() for k in KINDS}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            t = rec.info.get("SVTYPE")
            if t not in SVTYPES:
                continue
            fmt = rec.samples[0]
            gt = fmt.get("GT")
            if gt is None or not any(a == 1 for a in gt if a is not None):
                continue
            flt = set(rec.filter.keys())
            if flt == {"PASS"} or not flt:
                kind = "pass"
            elif flt == {"LowSQ"} and float(fmt.get("SQ") or 0.0) >= rescue_min_sq:
                kind = "rescuable"
            else:
                continue
            sc = size_class(rec.stop - rec.start)
            out[kind][(t, sc)] += 1
            out[kind][(t, "all")] += 1
    return out


@dataclass
class PreStats:
    name: str
    counts: dict[str, dict[str, Counter]]  # sample -> kind -> Counter
    wave: dict[str, float] = field(default_factory=dict)

    @property
    def samples(self) -> list[str]:
        return sorted(self.counts)

    def values(self, kind: str, key: tuple[str, str], only: set[str] | None = None) -> list[int]:
        return [c[kind][key] for s, c in sorted(self.counts.items()) if only is None or s in only]


def run_samples(run_dir: str | Path) -> list[tuple[str, Path]]:
    d = Path(run_dir) / "samples"
    out = []
    for p in sorted(d.iterdir()) if d.is_dir() else []:
        v = p / f"{p.name}.filtered.vcf.gz"
        if v.exists():
            out.append((p.name, v))
    return out


def _count_task(task):
    sample, path, rescue_min_sq = task
    return sample, count_sample(path, rescue_min_sq)


def load_filtered(
    run_dir: str | Path, name: str | None = None, rescue_min_sq: float = 1.0, jobs: int = 1,
    wave_from: str | Path | None = None,
) -> PreStats:
    """Counts for every ``samples/<S>/<S>.filtered.vcf.gz`` of the run; wave
    scores from its stores (or ``wave_from``'s)."""
    from ..parallel import pmap

    found = run_samples(run_dir)
    if not found:
        raise FileNotFoundError(f"no samples/<S>/<S>.filtered.vcf.gz under {run_dir}")
    counts = dict(pmap(_count_task, [(s, p, rescue_min_sq) for s, p in found], jobs))
    wave = read_waves(wave_from or run_dir, list(counts))
    return PreStats(name=name or Path(run_dir).name, counts=counts, wave=wave)


def pregate(
    baseline: PreStats,
    candidate: PreStats,
    max_increase: float = 0.25,
    max_increase_rescuable: float = 0.5,
    min_abs: float = 0.5,
) -> tuple[list[dict], list[str]]:
    """Rows and alarms: mean per genome of PASS and rescuable calls, over the
    genomes common to both runs, then within each wave stratum of >=
    MIN_STRATUM of them. An alarm needs a relative rise > max_increase
    (max_increase_rescuable for the rescuable calls) AND an absolute rise >=
    min_abs calls per genome.

    The rescuable limit is looser: those calls only reach the cohort when they
    recur, and a deployable candidate (cohort gate OK in every wave stratum)
    lifted the wavy genomes' rescuable DUP 1-10 kb by 30-36 %, while every
    candidate the cohort gate rejected moved a count by >= 52 % (most by 100 %
    to 1,000 %)."""
    common = set(baseline.counts) & set(candidate.counts)
    waves = {**baseline.wave, **candidate.wave}
    groups: list[tuple[str, set[str]]] = [("all", common)]
    groups += [(n, m & common) for n, m in wave_strata(waves) if len(m & common) >= MIN_STRATUM]
    rows, alarms = [], []
    for wname, only in groups:
        for kind in KINDS:
            for key in KEYS:
                bv, cv = baseline.values(kind, key, only), candidate.values(kind, key, only)
                b = statistics.fmean(bv) if bv else 0.0
                c = statistics.fmean(cv) if cv else 0.0
                lim = max_increase if kind == "pass" else max_increase_rescuable
                flag = _rise(b, c, lim, min_abs)
                if flag:
                    where = "" if wname == "all" else f" [wave {wname}, {len(only)} genomes]"
                    alarms.append(f"{kind} {key[0]} {key[1]}{where}: mean per genome {b:.2f} -> {c:.2f}")
                rows.append({
                    "run": candidate.name, "wave": wname, "n_genomes": len(only), "kind": kind,
                    "svtype": key[0], "size": key[1], "mean_per_genome": round(c, 3),
                    "median_per_genome": statistics.median(cv) if cv else 0, "total": sum(cv),
                    "baseline_mean": round(b, 3), "change": f"{c - b:+.2f}",
                    "flag": "ALARM" if flag else "",
                })
    return rows, alarms


def summary_lines(s: PreStats) -> list[str]:
    out = [f"== {s.name}: {len(s.counts)} genomes" + (
        f" (wave score missing for {len(s.counts) - len(s.wave)})" if len(s.wave) < len(s.counts) else "")]
    for kind in KINDS:
        parts = []
        for t in SVTYPES:
            v = {sc: sum(s.values(kind, (t, sc))) for sc in SIZE_CLASSES}
            parts.append(f"{t} " + " / ".join(f"{v[sc]}" for sc in SIZE_CLASSES))
        out.append(f"   {kind:>9} totals (<1kb / 1-10kb / >=10kb): " + "; ".join(parts))
    return out


def write_table(rows: list[dict], out: str | Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
