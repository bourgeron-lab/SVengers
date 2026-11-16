"""Cohort site correspondence between two runs (issue #9)."""

import csv

import pysam

from svengers.cohort.sitemap import map_sites, write_map


def _cohort_vcf(path, sites, spec=True):
    """(id, contig, start0, end, svtype) sites; spec or legacy (pre-1.4.0) coordinates."""
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=svengers")
    if spec:
        hdr.add_line("##svengers_coords=padded")
    hdr.contigs.add("chr1", length=10_000_000)
    for k, d in (("DEL", "Deletion"), ("DUP", "Duplication")):
        hdr.add_line(f'##ALT=<ID={k},Description="{d}">')
    hdr.info.add("SVTYPE", "1", "String", "type")
    hdr.info.add("END", "1", "Integer", "end")
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for sid, contig, s, e, t in sites:
            # spec: POS = padding base (start0), END = last base; legacy: both one to the right
            pos0, stop = (s - 1, e) if spec else (s, e + 1)
            rec = vf.new_record(contig=contig, start=pos0, stop=stop, alleles=("N", f"<{t}>"), id=sid)
            rec.info["SVTYPE"] = t
            vf.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_site_map_relations(tmp_path):
    """A validation keyed on a site follows the event to the next merge: the
    same site, a moved one, one split in two, two merged in one, a lost and a
    new one; legacy and spec coordinates give the same keys."""
    old = tmp_path / "old.vcf.gz"
    new = tmp_path / "new.vcf.gz"
    _cohort_vcf(old, [
        ("o1", "chr1", 10_000, 20_000, "DEL"),      # same
        ("o2", "chr1", 100_000, 110_000, "DEL"),    # moved
        ("o3", "chr1", 200_000, 300_000, "DUP"),    # split
        ("o4", "chr1", 400_000, 450_000, "DEL"),    # merged with o5
        ("o5", "chr1", 401_000, 452_000, "DEL"),
        ("o6", "chr1", 600_000, 610_000, "DEL"),    # lost
        ("o7", "chr1", 700_000, 710_000, "DUP"),    # a DEL takes its place: lost
    ], spec=False)
    _cohort_vcf(new, [
        ("n1", "chr1", 10_000, 20_000, "DEL"),
        ("n2", "chr1", 100_500, 110_200, "DEL"),
        ("n3", "chr1", 200_000, 260_000, "DUP"),
        ("n4", "chr1", 205_000, 300_000, "DUP"),
        ("n5", "chr1", 400_500, 451_000, "DEL"),
        ("n6", "chr1", 700_000, 710_000, "DEL"),    # new
    ])
    rows, counts = map_sites(old, new)
    assert counts == {"same": 1, "moved": 1, "split": 2, "merged": 2, "lost": 2, "new": 1}
    rel = {(r.get("old_site_id"), r.get("new_site_id")): r["relation"] for r in rows}
    assert rel[("o1", "n1")] == "same" and rel[("o2", "n2")] == "moved"
    assert rel[("o3", "n3")] == rel[("o3", "n4")] == "split"
    assert rel[("o4", "n5")] == rel[("o5", "n5")] == "merged"
    assert rel[("o6", None)] == "lost" and rel[(None, "n6")] == "new"
    same = next(r for r in rows if r.get("old_site_id") == "o1")
    assert same["old_site_key"] == same["new_site_key"] == "chr1:10001-20000:del"

    out = tmp_path / "map.tsv"
    write_map(rows, out)
    got = list(csv.DictReader(open(out), delimiter="\t"))
    assert len(got) == len(rows) and got[0]["old_site_id"] == "o1"
