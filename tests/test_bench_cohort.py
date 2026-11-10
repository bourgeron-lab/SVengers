"""Cohort gate: de novo burden per child against a baseline run."""

import csv

import pysam
from typer.testing import CliRunner

from vedette.bench.cohort import gate, load_run, size_class
from vedette.cli import app

REPORT = [
    "site_id", "contig", "start", "end", "svtype", "family_id", "child", "child_gt",
    "child_gq", "father_gt", "father_gq", "mother_gt", "mother_gq", "inheritance",
    "violation", "dnq",
]
CHILDREN = [f"K{i}" for i in range(6)]


def _run(tmp_path, name, sites, denovo, children=CHILDREN):
    """sites: [(id, svtype, start, end, filter)]; denovo: [(site_id, child, gqs)]."""
    d = tmp_path / name / "cohort"
    d.mkdir(parents=True)
    h = pysam.VariantHeader()
    h.contigs.add("chr1", length=10_000_000)
    h.info.add("SVTYPE", "1", "String", "t")
    h.filters.add("LowQual", None, None, "low")
    with pysam.VariantFile(str(d / "cohort.vcf.gz"), "wz", header=h) as out:
        for sid, t, s, e, flt in sorted(sites, key=lambda x: x[2]):
            r = out.new_record(contig="chr1", start=s, stop=e, alleles=("N", f"<{t}>"), id=sid)
            r.info["SVTYPE"] = t
            r.filter.add(flt)
            out.write(r)
    by_id = {s[0]: s for s in sites}
    with open(d / "mendel_report.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT, delimiter="\t")
        w.writeheader()
        for sid, child, (gc, gf, gm) in denovo:
            _, t, s, e, _ = by_id[sid]
            w.writerow({
                "site_id": sid, "contig": "chr1", "start": s, "end": e, "svtype": t,
                "family_id": "F", "child": child, "child_gt": "0/1", "child_gq": gc,
                "father_gt": "0/0", "father_gq": gf, "mother_gt": "0/0", "mother_gq": gm,
                "inheritance": "denovo", "violation": 1, "dnq": "",
            })
            w.writerow({**{k: "" for k in REPORT}, "site_id": sid, "start": s, "end": e,
                        "svtype": t, "child": child, "inheritance": "biparental",
                        "child_gq": 50, "father_gq": 50, "mother_gq": 50})
    with open(d / "mendel_summary.tsv", "w") as fh:
        fh.write("family_id\tchild\tn_genotyped\tn_child_carrier\tn_denovo\t"
                 "n_violation\tviolation_rate\ttransmission_rate\n")
        for c in children:
            fh.write(f"F\t{c}\t10\t5\t1\t1\t0.02\t0.89\n")
    return tmp_path / name


SITES = [
    ("d_small", "DEL", 1_000, 1_500, "PASS"),
    ("d_mid", "DEL", 20_000, 25_000, "PASS"),
    ("d_low", "DEL", 40_000, 45_000, "LowQual"),
]


def _baseline(tmp_path):
    dn = [("d_small", c, (30, 30, 30)) for c in CHILDREN[:4]]
    dn += [("d_mid", "K0", (30, 10, 30))]             # parent GQ below the floor
    dn += [("d_low", c, (30, 30, 30)) for c in CHILDREN]  # non-PASS site
    return _run(tmp_path, "base", SITES, dn)


def test_size_class():
    assert size_class(999) == "<1kb"
    assert size_class(1_000) == "1-10kb"
    assert size_class(10_000) == ">=10kb"


def test_load_run_counts_pass_gq_denovo_per_child(tmp_path):
    s = load_run(_baseline(tmp_path))
    assert s.name == "base" and s.n_children == 6
    assert s.per_child[("all", "all")] == [1, 1, 1, 1, 0, 0]
    assert s.per_child[("DEL", "1-10kb")] == [0] * 6
    assert s.recurrent_sites == 0          # d_small: 4 children < 5
    assert s.pass_sites == {"DEL": 2}
    assert s.mvr == 0.02


def test_gate_flags_a_family_shared_dup_burst(tmp_path):
    base = load_run(_baseline(tmp_path))
    dups = [(f"u{i}", "DUP", 100_000 + i * 20_000, 112_000 + i * 20_000, "PASS")
            for i in range(4)]
    dn = [("d_small", c, (30, 30, 30)) for c in CHILDREN[:4]]
    dn += [(u[0], c, (40, 40, 40)) for u in dups for c in CHILDREN]
    cand = load_run(_run(tmp_path, "cand", SITES + dups, dn))
    rows, alarms = gate(base, cand)
    assert cand.median(("DUP", ">=10kb")) == 4 and cand.recurrent_sites == 4
    assert any("DUP >=10kb" in a for a in alarms)
    flagged = {(r["svtype"], r["size"]) for r in rows if r["flag"]}
    assert ("DEL", "<1kb") not in flagged
    # a candidate identical to the baseline passes
    assert gate(base, base)[1] == []


def test_gate_ignores_small_absolute_rises(tmp_path):
    base = load_run(_baseline(tmp_path))
    dn = [("d_small", c, (30, 30, 30)) for c in CHILDREN]  # median 1 -> 1, 0/1 in DEL <1kb
    cand = load_run(_run(tmp_path, "cand", SITES, dn))
    assert gate(base, cand)[1] == []


def test_wave_strata_catch_a_burst_in_the_wavy_children(tmp_path):
    """A burst confined to a minority of wavy children leaves the cohort
    median flat; the wave stratum raises the alarm."""
    import zarr

    from vedette.bench import cohort as bc

    kids = [f"K{i}" for i in range(30)]
    wavy = kids[:8]
    base = _run(tmp_path, "base", SITES, [("d_small", c, (30, 30, 30)) for c in kids[:4]], kids)
    dups = [(f"u{i}", "DUP", 100_000 + i * 20_000, 103_000 + i * 20_000, "PASS") for i in range(6)]
    dn = [(u[0], c, (40, 40, 40)) for u in dups for c in wavy]
    cand = _run(tmp_path, "cand", SITES + dups, dn, kids)
    for c in kids:  # only the baseline's stores carry the wave scores
        g = zarr.open_group(str(base / "samples" / c / f"{c}.evidence.zarr"), mode="w")
        g.attrs["wave_score"] = 0.25 if c in wavy else 0.09
    old = bc.MIN_STRATUM
    bc.MIN_STRATUM = 5
    try:
        b, k = load_run(base), load_run(cand)
        assert b.wave["K0"] == 0.25 and k.wave == {}
        assert k.median(("DUP", "1-10kb")) == 0          # 8 of 30 children
        rows, alarms = gate(b, k)
        dup = [a for a in alarms if "DUP" in a]
        assert dup and all("[wave >=0.20, 8 children]" in a for a in dup)
        assert not any(r["flag"] for r in rows
                       if r["wave"] in ("all", "<0.128") and r["svtype"] != "recurrent_sites")
        assert any(">=0.20 (8 children)" in line for line in bc.summary_lines(b))
    finally:
        bc.MIN_STRATUM = old


def test_cli_exit_code(tmp_path):
    base = _baseline(tmp_path)
    dups = [(f"u{i}", "DUP", 100_000 + i * 20_000, 112_000 + i * 20_000, "PASS")
            for i in range(4)]
    dn = [(u[0], c, (40, 40, 40)) for u in dups for c in CHILDREN]
    cand = _run(tmp_path, "cand", SITES + dups, dn)
    r = CliRunner().invoke(app, ["bench", "cohort", str(base), f"new={cand}",
                                 "--out", str(tmp_path / "t.tsv"), "--jobs", "1"])
    assert r.exit_code == 1, r.output
    assert "new vs base: ALARM" in r.output
    r = CliRunner().invoke(app, ["bench", "cohort", str(base), str(base), "--jobs", "1"])
    assert r.exit_code == 0, r.output
