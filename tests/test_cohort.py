import numpy as np
import pytest

from vedette.cohort.merge import MemberCall, SiteCluster, _compatible, cluster_calls
from vedette.cohort.metrics import hwe_exact


def _mc(sample, start, end, precise=True, gq=99):
    return MemberCall(sample=sample, start=start, end=end, precise=precise, gq=gq)


def test_identical_calls_cluster():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 1010, 5005), "DEL")]
    clusters = cluster_calls(calls)
    assert len(clusters) == 1
    assert clusters[0].members[0].sample == "A"
    assert len(clusters[0].members) == 2


def test_del_dup_never_merge():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 1000, 5000), "DUP")]
    assert len(cluster_calls(calls)) == 2


def test_distant_calls_stay_separate():
    calls = [(_mc("A", 1000, 5000), "DEL"), (_mc("B", 50_000, 55_000), "DEL")]
    assert len(cluster_calls(calls)) == 2


def test_no_transitive_chaining():
    """A drifting chain must not fuse into one mega-cluster (AF inflation guard)."""
    calls = []
    start, end = 10_000, 20_000
    for i in range(30):
        calls.append((_mc(f"S{i}", start, end, precise=False), "DEL"))
        start += 1_500  # each next call shifts; pairwise-adjacent are compatible
        end += 1_500
    clusters = cluster_calls(sorted(calls, key=lambda t: t[0].start))
    assert len(clusters) >= 5  # centroid rule caps the diameter
    assert max(len(c.members) for c in clusters) < 15


def test_precise_needs_tight_breakpoints():
    cluster = SiteCluster(svtype="DEL", members=[_mc("A", 10_000, 60_000)])
    assert _compatible(_mc("B", 10_050, 60_100), cluster)          # within tolerance
    assert not _compatible(_mc("C", 14_000, 64_000), cluster)      # RO ok, bp too far


def test_hwe_extremes():
    # perfect equilibrium-ish: common variant
    assert hwe_exact(50, 25, 25) > 0.01
    # gross excess of hets (typical artifact signature)
    assert hwe_exact(100, 0, 0) < 1e-6
    # monomorphic
    assert hwe_exact(0, 0, 100) == 1.0


def _hwe_brute_force(het, hom_alt, hom_ref):
    """Exact conditional test by direct enumeration (Wigginton definition)."""
    from math import exp, lgamma

    n = het + hom_alt + hom_ref
    rare = 2 * min(hom_alt, hom_ref) + het

    def log_prob(h):
        if (rare - h) % 2:
            return None
        r = (rare - h) // 2
        c = n - h - r
        if r < 0 or c < 0:
            return None
        # P(h | n, rare) ∝ n! / (r! h! c!) * 2^h
        return (
            lgamma(n + 1) - lgamma(r + 1) - lgamma(h + 1) - lgamma(c + 1)
            + h * np.log(2)
        )

    logs = {h: lp for h in range(rare + 1) if (lp := log_prob(h)) is not None}
    mx = max(logs.values())
    probs = {h: exp(lp - mx) for h, lp in logs.items()}
    total = sum(probs.values())
    obs = probs[het]
    return sum(p for p in probs.values() if p <= obs * 1.0000001) / total


def test_hwe_matches_brute_force():
    for het, hom_alt, hom_ref in [(5, 2, 13), (57, 13, 298), (10, 0, 90), (3, 3, 3)]:
        assert hwe_exact(het, hom_alt, hom_ref) == pytest.approx(
            _hwe_brute_force(het, hom_alt, hom_ref), rel=1e-6
        )


def _write_filtered_vcf(path, sample, records):
    """records: (contig, start, end, svtype, filters, sq)"""
    import pysam

    header = pysam.VariantHeader()
    header.contigs.add("chr22", length=50_818_468)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("END", "1", "Integer", "End")
    header.filters.add("LowSQ", None, None, "below SQ threshold")
    header.filters.add("MaxN", None, None, "N overlap")
    header.formats.add("GT", "1", "String", "Genotype")
    header.formats.add("GQ", "1", "Integer", "Genotype quality")
    header.formats.add("SQ", "1", "Float", "Calibrated site quality")
    header.add_sample(sample)
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, end, svtype, filters, sq in sorted(records):
            rec = out.new_record(
                contig=contig, start=start, stop=end,
                alleles=("N", f"<{svtype}>"),
            )
            rec.info["SVTYPE"] = svtype
            for f in filters:
                rec.filter.add(f)
            rec.samples[sample]["GT"] = (0, 1)
            rec.samples[sample]["GQ"] = 99
            rec.samples[sample]["SQ"] = sq
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_merge_rescues_recurrent_subthreshold_sites(tmp_path):
    import pysam

    from vedette.cohort.merge import run_merge
    from vedette.reference import Contig

    # 5 carriers of the same LowSQ 30 kb DEL (the SHANK3 family case),
    # 1 singleton LowSQ elsewhere (must NOT be rescued),
    # 1 PASS call elsewhere (control)
    paths = []
    for i in range(5):
        p = tmp_path / f"s{i}.vcf.gz"
        recs = [("chr22", 50_654_000, 50_684_001, "DEL", ["LowSQ"], 2.1)]
        if i == 0:
            recs.append(("chr22", 10_000_000, 10_030_000, "DEL", ["LowSQ"], 3.0))
            recs.append(("chr22", 20_000_000, 20_030_000, "DEL", [], 20.0))
            recs.append(("chr22", 30_000_000, 30_030_000, "DEL", ["MaxN"], 15.0))
        _write_filtered_vcf(p, f"s{i}", recs)
        paths.append(p)

    out = tmp_path / "sites.vcf.gz"
    contigs = [Contig("chr22", 50_818_468, "autosome")]
    n = run_merge(paths, contigs, out)
    # read-back stop is END = POS + |SVLEN| (htslib symbolic-allele rule),
    # one base past the half-open end — the codebase-wide convention
    sites = {(r.start, r.stop - 1): r for r in pysam.VariantFile(str(out))}
    # rescued familial site present and flagged
    assert (50_654_000, 50_684_001) in sites
    assert sites[(50_654_000, 50_684_001)].info["RESCUE"] == "RECURRENT"
    assert sites[(50_654_000, 50_684_001)].info["NSAMPLES"] == 5
    # PASS control present without flag; singleton LowSQ and hard-filter absent
    assert (20_000_000, 20_030_000) in sites
    assert "RESCUE" not in sites[(20_000_000, 20_030_000)].info
    assert (10_000_000, 10_030_000) not in sites
    assert (30_000_000, 30_030_000) not in sites
    assert n == 2
    # per-sample best SQ survives the merge (SQMAX/SQMEAN)
    assert sites[(20_000_000, 20_030_000)].info["SQMAX"] == 20.0
    assert abs(sites[(50_654_000, 50_684_001)].info["SQMAX"] - 2.1) < 1e-6
    assert abs(sites[(50_654_000, 50_684_001)].info["SQMEAN"] - 2.1) < 1e-6


def _write_sites_vcf(path, records):
    """records: (site_id, contig, start, end, svtype, nsamples)"""
    import pysam

    header = pysam.VariantHeader()
    header.contigs.add("chr22", length=50_818_468)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("END", "1", "Integer", "End")
    header.info.add("NSAMPLES", "1", "Integer", "carriers")
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for sid, contig, start, end, svtype, ns in records:
            rec = out.new_record(contig=contig, start=start, stop=end,
                                 alleles=("N", f"<{svtype}>"), id=sid)
            rec.info["SVTYPE"] = svtype
            rec.info["NSAMPLES"] = ns
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def _write_reference_cohort(path, records):
    """records: (contig, start, end, svtype, af, ns, ncarrier, filters)"""
    import pysam

    header = pysam.VariantHeader()
    header.contigs.add("chr22", length=50_818_468)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("END", "1", "Integer", "End")
    header.info.add("AF", "A", "Float", "AF")
    header.info.add("NS", "1", "Integer", "NS")
    header.info.add("NCARRIER", "1", "Integer", "carriers")
    header.filters.add("CohortArtifact", None, None, "artifact")
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, end, svtype, af, ns, ncar, filters in records:
            rec = out.new_record(contig=contig, start=start, stop=end,
                                 alleles=("N", f"<{svtype}>"))
            rec.info["SVTYPE"] = svtype
            rec.info["AF"] = (af,)
            rec.info["NS"] = ns
            rec.info["NCARRIER"] = ncar
            for f in filters:
                rec.filter.add(f)
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_gtr_discordance_is_directional(tmp_path):
    """Re-genotyping finding MORE carriers than the callers nominated is the
    point of genotyping from the stores (it recovered the causal SHANK3
    deletion on the PMS cohort: 6 nominators, 64 carriers). Only the other
    direction — nominations that re-genotyping fails to reproduce — is a
    red flag."""
    import csv

    import pysam

    from vedette.cohort.genotyper import GENOTYPE_COLUMNS
    from vedette.cohort.metrics import run_cohort_filter

    # UP: 3 nominators, 6 of 8 samples carry after re-genotyping (sensitivity)
    # DOWN: 8 nominators, only 1 carrier reproduced (irreproducible calls)
    sites = [("UP", "chr22", 1_000_000, 1_030_000, "DEL", 3),
             ("DOWN", "chr22", 2_000_000, 2_030_000, "DEL", 8)]
    _write_sites_vcf(tmp_path / "sites.vcf.gz", sites)
    gt_paths = []
    for k in range(8):
        rows = []
        for sid, contig, start, end, svtype, _ in sites:
            if sid == "UP":
                gt, cn = ("0/1", 1) if k < 6 else ("0/0", 2)
            else:
                gt, cn = ("0/1", 1) if k == 0 else ("0/0", 2)
            rows.append(dict(zip(GENOTYPE_COLUMNS, [
                sid, contig, start, end, svtype, f"S{k}", gt, 40, cn, 0.99, "", 0, 0, 0, 0, "",
            ])))
        path = tmp_path / f"S{k}.gt.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        gt_paths.append(path)

    run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "cohort.vcf.gz")
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "cohort.vcf.gz"))}
    assert "CohortArtifact" not in out["UP"].filter.keys()
    assert "CohortArtifact" in out["DOWN"].filter.keys()
    # the absolute discordance is still reported for both
    assert out["UP"].info["GTR_DISC"] == pytest.approx(1.0)


def test_cohort_filter_post_processing_flags(tmp_path):
    import csv

    import pysam

    from vedette.cohort.genotyper import GENOTYPE_COLUMNS
    from vedette.cohort.metrics import run_cohort_filter

    # A: DEL carried 1/1 by all 6 samples (reference artifact, HWE-blind)
    # B: DUP carried by 4/6, every carrier at CN=1 (type/depth discordant)
    # C: DEL carried by 1/6, matches a reference-cohort site carried by 92 %
    # D: DEL carried by 2/6 with consistent CN, unknown to the reference (clean)
    sites = [("A", "chr22", 1_000_000, 1_001_000, "DEL", 6),
             ("B", "chr22", 2_000_000, 2_000_300, "DUP", 4),
             ("C", "chr22", 3_000_000, 3_010_000, "DEL", 1),
             ("D", "chr22", 4_000_000, 4_020_000, "DEL", 2)]
    _write_sites_vcf(tmp_path / "sites.vcf.gz", sites)
    _write_reference_cohort(tmp_path / "ref.vcf.gz", [
        ("chr22", 3_000_050, 3_009_900, "DEL", 0.48, 2000, 1840, ["CohortArtifact"]),
        ("chr22", 4_000_000, 4_020_000, "DUP", 0.9, 2000, 1900, []),  # wrong type: no match
    ])
    gt_paths = []
    for k in range(6):
        smp = f"S{k}"
        rows = []
        for sid, contig, start, end, svtype, _ in sites:
            if sid == "A":
                gt, cn = "1/1", 0
            elif sid == "B":
                gt, cn = ("0/1", 1) if k < 4 else ("0/0", 2)
            elif sid == "C":
                gt, cn = ("0/1", 1) if k == 0 else ("0/0", 2)
            else:
                gt, cn = ("0/1", 1) if k < 2 else ("0/0", 2)
            rows.append(dict(zip(GENOTYPE_COLUMNS, [
                sid, contig, start, end, svtype, smp, gt, 40, cn, 0.99, "", 0, 0, 0, 0, "",
            ])))
        path = tmp_path / f"{smp}.gt.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        gt_paths.append(path)

    counts = run_cohort_filter(
        tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "cohort.vcf.gz",
        known_sites=tmp_path / "ref.vcf.gz",
    )
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "cohort.vcf.gz"))}
    assert "HighCarrierFrac" in out["A"].filter.keys()
    assert out["A"].info["CARRIER_FRAC"] == 1.0
    assert "CNDiscordant" in out["B"].filter.keys()
    assert out["B"].info["CN_DISC"] == 1.0
    assert "HighCarrierFrac" not in out["B"].filter.keys()
    assert "CommonRef" in out["C"].filter.keys()
    assert abs(out["C"].info["REF_CARRIER_FRAC"] - 0.92) < 1e-6
    assert out["C"].info["REF_FILTER"] == "CohortArtifact"
    assert list(out["D"].filter.keys()) == ["PASS"]
    assert "REF_AF" not in out["D"].info
    assert counts["high_carrier"] == 1 and counts["cn_discordant"] == 1
    assert counts["common_ref"] == 1 and counts["pass"] == 1
