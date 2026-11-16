"""Pre-gate: PASS and rescuable calls per genome against a baseline, on filter-only runs."""

from collections import Counter

import pysam
from typer.testing import CliRunner

from svengers.bench.pregate import PreStats, count_sample, load_filtered, pregate
from svengers.cli import app


def _vcf(path, records):
    """records: [(svtype, start, end, filter, gt, sq)]."""
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=10_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.info.add("END", "1", "Integer", "e")
    h.filters.add("LowSQ", None, None, "low SQ")
    h.filters.add("NoEvidence", None, None, "no evidence")
    h.formats.add("GT", "1", "String", "genotype")
    h.formats.add("SQ", "1", "Float", "score")
    h.add_sample("S")
    with pysam.VariantFile(str(path), "wz", header=h) as out:
        for t, s, e, flt, gt, sq in sorted(records, key=lambda r: r[1]):
            r = out.new_record(contig="chr1", start=s, stop=e, alleles=("N", f"<{t}>"))
            r.info["SVTYPE"] = t
            for f in flt.split(";"):
                r.filter.add(f)
            r.samples["S"]["GT"] = gt
            r.samples["S"]["SQ"] = sq
            out.write(r)


def _run(tmp_path, name, per_sample):
    """per_sample: {sample: records}; Nextflow layout samples/<S>/<S>.filtered.vcf.gz."""
    for s, recs in per_sample.items():
        d = tmp_path / name / "samples" / s
        d.mkdir(parents=True)
        _vcf(d / f"{s}.filtered.vcf.gz", recs)
    return tmp_path / name


BASE = [
    ("DEL", 1_000, 1_400, "PASS", (0, 1), 9.0),
    ("DEL", 20_000, 25_000, "PASS", (1, 1), 8.0),
    ("DUP", 40_000, 42_000, "LowSQ", (0, 1), 0.4),
]


def test_count_sample_pass_and_rescuable(tmp_path):
    p = tmp_path / "s.vcf.gz"
    _vcf(p, BASE + [
        ("DUP", 50_000, 52_000, "LowSQ", (0, 1), 1.5),              # rescuable
        ("DUP", 60_000, 75_000, "LowSQ", (0, 1), 1.0),              # rescuable, at the floor
        ("DUP", 80_000, 82_000, "NoEvidence;LowSQ", (0, 1), 3.0),   # hard filter: never rescued
        ("DUP", 90_000, 92_000, "PASS", (0, 0), 9.0),               # non-carrier
    ])
    c = count_sample(p)
    assert c["pass"] == Counter({("DEL", "<1kb"): 1, ("DEL", "1-10kb"): 1, ("DEL", "all"): 2})
    assert c["rescuable"] == Counter({("DUP", "1-10kb"): 1, ("DUP", ">=10kb"): 1, ("DUP", "all"): 2})
    assert count_sample(p, rescue_min_sq=2.0)["rescuable"] == Counter()


def test_pregate_flags_a_dup_flood_not_a_stable_run(tmp_path):
    samples = [f"S{i}" for i in range(4)]
    base = load_filtered(_run(tmp_path, "base", {s: BASE for s in samples}))
    same = load_filtered(_run(tmp_path, "same", {s: BASE for s in samples}))
    flood = BASE + [("DUP", 100_000 + 10_000 * k, 101_500 + 10_000 * k, "PASS", (0, 1), 5.0) for k in range(2)]
    cand = load_filtered(_run(tmp_path, "cand", {s: flood for s in samples}))
    assert len(base.counts) == 4 and base.name == "base"
    _, alarms = pregate(base, same)
    assert alarms == []
    rows, alarms = pregate(base, cand)
    assert alarms == [
        "pass DUP all: mean per genome 0.00 -> 2.00",
        "pass DUP 1-10kb: mean per genome 0.00 -> 2.00",
    ]
    r = next(r for r in rows if r["kind"] == "pass" and r["svtype"] == "DUP" and r["size"] == "1-10kb")
    assert (r["total"], r["baseline_mean"], r["flag"]) == (8, 0.0, "ALARM")


def _stats(name, n, rescuable, wave):
    counts = {}
    for i in range(n):
        k = rescuable(i)
        counts[f"S{i}"] = {"pass": Counter(), "rescuable": Counter({("DUP", "1-10kb"): k, ("DUP", "all"): k})}
    return PreStats(name=name, counts=counts, wave={f"S{i}": wave(i) for i in range(n)})


def test_rescuable_limit_is_looser_and_strata_catch_wavy_genomes():
    wave = lambda i: 0.25 if i < 10 else 0.05  # noqa: E731 - 10 wavy genomes of 40
    base = _stats("base", 40, lambda i: 10, wave)
    # +40 % everywhere: under the rescuable limit (50 %), over the PASS one
    _, alarms = pregate(base, _stats("mild", 40, lambda i: 14, wave))
    assert alarms == []
    _, alarms = pregate(base, _stats("mild", 40, lambda i: 14, wave), max_increase_rescuable=0.25)
    assert "rescuable DUP 1-10kb: mean per genome 10.00 -> 14.00" in alarms
    # x3 in the wavy genomes only: +50 % on the mean of 40, flagged within the stratum
    _, alarms = pregate(base, _stats("wavy", 40, lambda i: 30 if i < 10 else 10, wave))
    assert "rescuable DUP 1-10kb [wave >=0.20, 10 genomes]: mean per genome 10.00 -> 30.00" in alarms
    assert "rescuable DUP 1-10kb: mean per genome 10.00 -> 15.00" not in alarms


def test_cli_exit_code(tmp_path):
    samples = ["A", "B"]
    base = _run(tmp_path, "base", {s: BASE for s in samples})
    flood = BASE + [("DUP", 100_000 + 10_000 * k, 101_500 + 10_000 * k, "PASS", (0, 1), 5.0) for k in range(3)]
    cand = _run(tmp_path, "cand", {s: flood for s in samples})
    out = tmp_path / "pre.tsv"
    r = CliRunner().invoke(app, ["bench", "pregate", f"prod={base}", f"new={base}", "--jobs", "1"])
    assert r.exit_code == 0, r.output
    assert "-- new vs prod: OK" in r.output
    r = CliRunner().invoke(app, ["bench", "pregate", f"prod={base}", f"new={cand}", "--out", str(out), "--jobs", "1"])
    assert r.exit_code == 1
    assert "-- new vs prod: ALARM" in r.output and "pass DUP 1-10kb: mean per genome 0.00 -> 3.00" in r.output
    assert out.read_text().splitlines()[0].startswith("run\twave\tn_genomes\tkind")
