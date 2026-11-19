import numpy as np
import pytest

from svengers.cohort.merge import MemberCall, SiteCluster, _compatible, cluster_calls
from svengers.cohort.metrics import hwe_exact
from svengers.vcfout import add_coords_line, span, sv_record


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


def test_sweep_matches_the_exhaustive_scan():
    """cluster_calls drops clusters that can no longer take a member from its
    scan and keeps running medians: same clusters, same order as scanning every
    cluster within the 1 Mb horizon with medians recomputed (v2.1.0)."""
    def reference(calls):
        open_, done = {"DEL": [], "DUP": []}, []
        for call, t in calls:
            keep = []
            for members in open_[t]:
                (done if _centroid(members)[1] + 1_000_000 < call.start else keep).append(members)
            open_[t] = keep
            best, best_ro = None, 0.0
            for members in keep:
                cl = SiteCluster(svtype=t, members=list(members))
                if _compatible(call, cl):
                    s, e = _centroid(members)
                    ro = (min(call.end, e) - max(call.start, s)) / max(call.end - call.start, e - s)
                    if ro > best_ro:
                        best, best_ro = members, ro
            if best is None:
                keep.append([call])
            else:
                best.append(call)
        return done + open_["DEL"] + open_["DUP"]

    def _centroid(members):
        return (int(np.median([m.start for m in members])), int(np.median([m.end for m in members])))

    rng = np.random.default_rng(7)
    calls = []
    for i in range(1_500):
        start = int(rng.integers(0, 3_000_000))
        size = int(rng.choice([300, 2_000, 9_000, 80_000, 1_500_000]) * rng.uniform(0.9, 1.1))
        if i % 3 == 0 and calls:   # recurrent sites with jittered breakpoints
            m = calls[int(rng.integers(len(calls)))][0]
            start, size = m.start + int(rng.integers(-200, 200)), m.end - m.start + int(rng.integers(-200, 200))
        calls.append((_mc(f"S{i}", start, start + max(size, 60), precise=bool(i % 2)),
                       "DEL" if i % 5 else "DUP"))
    # a long cluster created before a short one, both closed after a > 1 Mb gap
    for name, start, end in (("Y", 5_000_000, 5_050_000), ("X", 5_001_000, 5_003_000),
                             ("Z", 5_010_000, 5_012_000), ("W", 7_000_000, 7_001_000)):
        calls.append((_mc(name, start, end), "DEL"))
    calls.sort(key=lambda t: t[0].start)
    got = [[(m.sample, m.start, m.end) for m in c.members] for c in cluster_calls(calls)]
    want = [[(m.sample, m.start, m.end) for m in c] for c in reference(calls)]
    assert got == want
    assert max(len(c) for c in got) > 5


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


def _write_filtered_vcf(path, sample, records, contigs=(("chr22", 50_818_468),)):
    """records: (contig, start, end, svtype, filters, sq)"""
    import pysam

    header = pysam.VariantHeader()
    for name, length in contigs:
        header.contigs.add(name, length=length)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    add_coords_line(header)
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("SVLEN", "1", "Integer", "Signed length")
    header.info.add("END", "1", "Integer", "End")
    header.filters.add("LowSQ", None, None, "below SQ threshold")
    header.filters.add("MaxN", None, None, "N overlap")
    header.formats.add("GT", "1", "String", "Genotype")
    header.formats.add("GQ", "1", "Integer", "Genotype quality")
    header.formats.add("SQ", "1", "Float", "Calibrated site quality")
    header.add_sample(sample)
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, end, svtype, filters, sq in sorted(records):
            rec = sv_record(out, contig, start, end, svtype)
            for f in filters:
                rec.filter.add(f)
            rec.samples[sample]["GT"] = (0, 1)
            rec.samples[sample]["GQ"] = 99
            rec.samples[sample]["SQ"] = sq
            out.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_merge_rescues_recurrent_subthreshold_sites(tmp_path):
    import pysam

    from svengers.cohort.merge import run_merge
    from svengers.reference import Contig

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
    sites = {span(r): r for r in pysam.VariantFile(str(out))}
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


def test_merge_clamps_site_end_to_contig(tmp_path):
    """Member calls from VCFs written before vcfout clamped END overhang the
    contig by up to one depth bin; the cohort site must still end on it."""
    import pysam

    from svengers.cohort.merge import run_merge
    from svengers.reference import Contig

    # the published dataset_paper per-sample VCFs: END one 1 kb bin past chr12
    chr12_len = 133_275_309
    contigs = (("chr12", chr12_len),)
    paths = []
    for i in range(3):
        p = tmp_path / f"s{i}.vcf.gz"
        recs = [("chr12", 133_250_000, 133_276_001, "DEL", [], 20.0)]
        _write_filtered_vcf(p, f"s{i}", recs, contigs=contigs)
        paths.append(p)

    out = tmp_path / "sites.vcf.gz"
    assert run_merge(paths, [Contig("chr12", chr12_len, "autosome")], out) == 1
    (site,) = pysam.VariantFile(str(out))
    assert span(site) == (133_250_000, chr12_len)
    assert site.stop == chr12_len  # was 133,276,002
    assert site.info["SVLEN"] == -(chr12_len - 133_250_000)


def _write_sites_vcf(path, records):
    """records: (site_id, contig, start, end, svtype, nsamples)"""
    import pysam

    header = pysam.VariantHeader()
    header.contigs.add("chr22", length=50_818_468)
    header.contigs.add("chrX", length=156_040_895)
    header.contigs.add("chrY", length=57_227_415)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    add_coords_line(header)
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("SVLEN", "1", "Integer", "Signed length")
    header.info.add("END", "1", "Integer", "End")
    header.info.add("NSAMPLES", "1", "Integer", "carriers")
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for sid, contig, start, end, svtype, ns in records:
            rec = sv_record(out, contig, start, end, svtype, sid)
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
    add_coords_line(header)
    header.info.add("SVTYPE", "1", "String", "SV type")
    header.info.add("SVLEN", "1", "Integer", "Signed length")
    header.info.add("END", "1", "Integer", "End")
    header.info.add("AF", "A", "Float", "AF")
    header.info.add("NS", "1", "Integer", "NS")
    header.info.add("NCARRIER", "1", "Integer", "carriers")
    header.filters.add("CohortArtifact", None, None, "artifact")
    header.filters.add("HighCarrierFrac", None, None, "carried by most samples")
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for contig, start, end, svtype, af, ns, ncar, filters in records:
            rec = sv_record(out, contig, start, end, svtype)
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

    from svengers.cohort.genotyper import GENOTYPE_COLUMNS
    from svengers.cohort.metrics import run_cohort_filter

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
            ], strict=False)))
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

    from svengers.cohort.genotyper import GENOTYPE_COLUMNS
    from svengers.cohort.metrics import run_cohort_filter

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
        # two FILTERs: they must not reach INFO joined by ';' (the INFO separator)
        ("chr22", 3_000_050, 3_009_900, "DEL", 0.48, 2000, 1840,
         ["CohortArtifact", "HighCarrierFrac"]),
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
            ], strict=False)))
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
    assert out["C"].info["REF_FILTER"] == "CohortArtifact|HighCarrierFrac"
    assert "HighCarrierFrac" not in out["C"].info
    assert list(out["D"].filter.keys()) == ["PASS"]
    assert "REF_AF" not in out["D"].info
    assert counts["high_carrier"] == 1 and counts["cn_discordant"] == 1
    assert counts["common_ref"] == 1 and counts["pass"] == 1


def test_cohort_filter_sex_chromosomes(tmp_path):
    """chrX: carrier fraction per sex, ploidy-relative CN discordance; the
    sexes come from the reference CN of each sample's chrX genotypes."""
    import csv

    import pysam

    from svengers.cohort.genotyper import GENOTYPE_COLUMNS
    from svengers.cohort.metrics import run_cohort_filter

    men, women = [f"M{k}" for k in range(6)], [f"F{k}" for k in range(6)]
    # X1: junction artifact, every woman 0/1 at CN 2 and half the men: 0.75 pooled
    # X2: DEL of three men (CN 0); X3: DUP of two men at CN 2
    sites = [("X1", "chrX", 10_000_000, 10_000_300, "DUP", 9),
             ("X2", "chrX", 20_000_000, 20_030_000, "DEL", 3),
             ("X3", "chrX", 30_000_000, 30_050_000, "DUP", 2)]
    sites += [(f"R{k}", "chrX", 40_000_000 + k * 100_000, 40_010_000 + k * 100_000, "DEL", 1)
              for k in range(12)]   # nobody carries them: the reference CN
    _write_sites_vcf(tmp_path / "sites.vcf.gz", sites)
    gt_paths = []
    for smp in men + women:
        male, k = smp[0] == "M", int(smp[1])
        ref = ("0/0", 1) if male else ("0/0", 2)
        rows = []
        for sid, contig, start, end, svtype, _ in sites:
            gt, cn = ref
            if sid == "X1" and (not male or k < 3):
                gt, cn = ("1/1", 1) if male else ("0/1", 2)
            elif sid == "X2" and male and k < 3:
                gt, cn = "1/1", 0
            elif sid == "X3" and male and k < 2:
                gt, cn = "1/1", 2
            rows.append(dict(zip(GENOTYPE_COLUMNS, [
                sid, contig, start, end, svtype, smp, gt, 40, cn, cn / 2, "", 0, 0, 0, 0, "",
            ], strict=False)))
        path = tmp_path / f"{smp}.gt.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        gt_paths.append(path)

    run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "cohort.vcf.gz")
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "cohort.vcf.gz"))}
    assert out["X1"].info["CARRIER_FRAC"] == 0.75
    assert "HighCarrierFrac" in out["X1"].filter.keys()        # 6/6 women
    assert list(out["X2"].filter.keys()) == ["PASS"]
    assert out["X3"].info["CN_DISC"] == 0.0                     # CN 2 is a gain for a man
    assert list(out["X3"].filter.keys()) == ["PASS"]


def test_cohort_filter_chry_no_depth(tmp_path):
    """chrY outside the PARs: a site whose carriers sit at their baseline copy
    number (junction-only), or shorter than one depth bin, is ChrYNoDepth; a
    deletion the CN supports and a pseudo-autosomal site are not."""
    import csv

    import pysam

    from svengers.cohort.genotyper import GENOTYPE_COLUMNS
    from svengers.cohort.metrics import run_cohort_filter

    men, women = [f"M{k}" for k in range(6)], [f"F{k}" for k in range(6)]
    # Y1: small DUP of four men at CN 1; Y2: DEL of two men at CN 0;
    # Y3: the same junction-only DUP inside PAR1; Y4: Y2 at 600 bp
    sites = [(f"R{k}", "chrX", 40_000_000 + k * 100_000, 40_010_000 + k * 100_000, "DEL", 1)
             for k in range(12)]   # nobody carries them: the reference CN
    sites += [("Y3", "chrY", 1_000_000, 1_000_300, "DUP", 4),
              ("Y1", "chrY", 10_000_000, 10_000_300, "DUP", 4),
              ("Y2", "chrY", 15_000_000, 15_005_000, "DEL", 2),
              ("Y4", "chrY", 16_000_000, 16_000_600, "DEL", 2)]
    _write_sites_vcf(tmp_path / "sites.vcf.gz", sites)
    gt_paths = []
    for smp in men + women:
        male, k = smp[0] == "M", int(smp[1])
        rows = []
        for sid, contig, start, end, svtype, _ in sites:
            if contig == "chrY":
                gt, cn = ("0/0", 1) if male else ("0/0", 0)
            else:
                gt, cn = ("0/0", 1) if male else ("0/0", 2)
            if sid in ("Y1", "Y3") and male and k < 4:
                gt, cn = "1/1", 1
            elif sid in ("Y2", "Y4") and male and k < 2:
                gt, cn = "1/1", 0
            rows.append(dict(zip(GENOTYPE_COLUMNS, [
                sid, contig, start, end, svtype, smp, gt, 40, cn, cn / 2, "", 0, 0, 0, 0, "",
            ], strict=False)))
        path = tmp_path / f"{smp}.gt.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        gt_paths.append(path)

    counts = run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "cohort.vcf.gz")
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "cohort.vcf.gz"))}
    assert list(out["Y1"].filter.keys()) == ["ChrYNoDepth"]
    assert list(out["Y2"].filter.keys()) == ["PASS"]
    assert list(out["Y3"].filter.keys()) == ["PASS"]
    assert list(out["Y4"].filter.keys()) == ["ChrYNoDepth"]   # CN 0, but 600 bp
    assert counts["chry_no_depth"] == 2


@pytest.mark.parametrize("child,father,mother,mode,male,expected", [
    # son: chrX from the mother only, coded hemizygous 1/1
    ((1, 1), (0, 0), (0, 1), "x", True, ("maternal", False)),
    ((1, 1), (1, 1), (0, 0), "x", True, ("denovo", True)),
    ((0, 0), (1, 1), (0, 1), "x", True, ("none", False)),
    ((0, 0), (0, 0), (1, 1), "x", True, ("none", True)),
    # daughter: always her father's chrX
    ((0, 1), (1, 1), (0, 0), "x", False, ("paternal", False)),
    ((0, 0), (1, 1), (0, 0), "x", False, ("none", True)),
    ((0, 1), (0, 0), (0, 0), "x", False, ("denovo", True)),
    # chrY: father to son
    ((1, 1), (1, 1), (0, 0), "y", True, ("paternal", False)),
    ((1, 1), (0, 0), (0, 0), "y", True, ("denovo", True)),
    # pseudo-autosomal: only "carrier child, both parents 0/0" is a violation
    ((1, 1), (0, 0), (0, 1), "par", True, ("maternal", False)),
    ((0, 0), (1, 1), (1, 1), "par", True, ("none", False)),
    # autosomes unchanged
    ((1, 1), (0, 0), (0, 1), "auto", True, ("maternal", True)),
])
def test_classify_trio_sex_chromosomes(child, father, mother, mode, male, expected):
    from svengers.cohort.mendel import classify_trio

    assert classify_trio(child, father, mother, mode, male) == expected


def test_site_mode():
    from svengers.cohort.mendel import site_mode

    assert site_mode("chr7", 0, 10) == "auto"
    assert site_mode("chrX", 1_000_000, 1_001_000) == "par"
    assert site_mode("chrX", 2_781_000, 2_790_000) == "par"     # overlaps PAR1
    assert site_mode("chrX", 50_000_000, 50_010_000) == "x"
    assert site_mode("chrY", 7_000_000, 7_010_000) == "y"


def test_cohort_filter_takes_the_reference_artifact_verdict(tmp_path):
    """A small cohort has neither the power nor the stability of the
    founders' HWE test: with --known-sites, a site common here whose
    reference site is CohortArtifact is CohortArtifact too (REF_ARTIFACT);
    a rare site keeps its own verdict."""
    import pysam

    from svengers.cohort.metrics import run_cohort_filter

    test_cohort_filter_post_processing_flags(tmp_path)    # C: AF 1/12, reference CohortArtifact
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "cohort.vcf.gz"))}
    assert "CohortArtifact" in out["C"].filter.keys() and out["C"].info.get("REF_ARTIFACT")
    assert "REF_ARTIFACT" not in out["D"].info
    gts = sorted(tmp_path.glob("S*.gt.tsv"))
    for kw in ({"artifact_af": 0.1}, {"known_artifacts": False}):    # C rare, or option off
        run_cohort_filter(tmp_path / "sites.vcf.gz", gts, tmp_path / "c2.vcf.gz",
                          known_sites=tmp_path / "ref.vcf.gz", **kw)
        c = {r.id: r for r in pysam.VariantFile(str(tmp_path / "c2.vcf.gz"))}["C"]
        assert "CohortArtifact" not in c.filter.keys() and "REF_ARTIFACT" not in c.info


def test_cohort_filter_trio_flags(tmp_path):
    import csv

    import pysam

    from svengers.cohort.genotyper import GENOTYPE_COLUMNS
    from svengers.cohort.metrics import run_cohort_filter
    from svengers.ped import parse_ped

    # 10 trios (P = father, M = mother, C = child)
    # G: fathers 0/1, half of the children carry      -> transmitted, PASS
    # T: fathers 0/1, no child carries                -> LowTransmission
    # N: 4 children carry, every parent 0/0, no depth -> DeNovoExcess
    # R: 2 children carry, parents 0/0 (NAHR hotspot) -> PASS
    # S: 4 children carry, every parent 0/0, RD 1.5   -> PASS (recurrent de novo with depth)
    fams = range(10)
    (tmp_path / "ped.tsv").write_text("".join(
        f"F{k}\tP{k}\t0\t0\t1\t1\nF{k}\tM{k}\t0\t0\t2\t1\nF{k}\tC{k}\tP{k}\tM{k}\t1\t2\n"
        for k in fams
    ))

    def carries(sid, smp):
        role, k = smp[0], int(smp[1:])
        if sid in ("G", "T"):
            return role == "P" or (sid == "G" and role == "C" and k < 5)
        return role == "C" and k < (2 if sid == "R" else 4)

    sites = [("G", "chr22", 1_000_000, 1_005_000, "DUP"),
             ("T", "chr22", 2_000_000, 2_005_000, "DUP"),
             ("N", "chr22", 3_000_000, 3_005_000, "DUP"),
             ("R", "chr22", 4_000_000, 4_005_000, "DUP"),
             ("S", "chr22", 5_000_000, 5_005_000, "DUP")]
    samples = [f"{r}{k}" for k in fams for r in "PMC"]
    _write_sites_vcf(tmp_path / "sites.vcf.gz", [
        (*s, sum(carries(s[0], x) for x in samples)) for s in sites
    ])
    gt_paths = []
    for smp in samples:
        path = tmp_path / f"{smp}.gt.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
            w.writeheader()
            for sid, contig, start, end, svtype in sites:
                c = carries(sid, smp)
                w.writerow(dict(zip(GENOTYPE_COLUMNS, [
                    sid, contig, start, end, svtype, smp, "0/1" if c else "0/0", 40,
                    3 if c else 2, 40, (1.05 if sid in "NR" else 1.5) if c else 1.0,
                    0, 0, "", 0, "",
                ], strict=False)))
        gt_paths.append(path)

    counts = run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "c.vcf.gz",
                               pedigree=parse_ped(tmp_path / "ped.tsv"))
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "c.vcf.gz"))}
    assert list(out["G"].filter.keys()) == ["PASS"]
    assert out["G"].info["TR_OPP"] == 10 and out["G"].info["TR_RATE"] == 0.5
    assert list(out["T"].filter.keys()) == ["LowTransmission"]
    assert out["T"].info["TR_RATE"] == 0.0
    assert list(out["N"].filter.keys()) == ["DeNovoExcess"]
    assert out["N"].info["DN_TRIO"] == 4 and out["N"].info["DN_FRAC"] == 1.0
    assert out["N"].info["DN_NODEPTH"] == 4
    assert list(out["R"].filter.keys()) == ["PASS"]
    assert list(out["S"].filter.keys()) == ["PASS"]
    assert out["S"].info["DN_TRIO"] == 4 and out["S"].info["DN_NODEPTH"] == 0
    assert counts["low_transmission"] == 1 and counts["denovo_excess"] == 1
    # FORMAT/FT: a carrier with no pair and RD 1.05 has no support of its own
    # (read as missing by the trio tools); RD 1.5 does; the site flags use every genotype
    from svengers.cohort.mendel import _gt_of
    assert out["N"].samples["C0"]["FT"] == "WeakSupport"
    assert _gt_of(out["N"], "C0") == (None, 0)
    assert out["G"].samples["C0"].get("FT") is None     # no weak genotype: no FT on the record
    assert _gt_of(out["G"], "C0") == ((0, 1), 40)
    assert out["N"].samples["C9"]["FT"] == "PASS"       # 0/0: never flagged
    assert out["N"].samples["P0"]["FT"] == "PASS"

    # thresholds: 11 trios needed -> T passes; 2 de novo children suffice -> R flagged
    run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "c2.vcf.gz",
                      pedigree=parse_ped(tmp_path / "ped.tsv"),
                      trio={"trio_min_opp": 11, "denovo_excess_n": 2})
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "c2.vcf.gz"))}
    assert list(out["T"].filter.keys()) == ["PASS"]
    assert list(out["R"].filter.keys()) == ["DeNovoExcess"]

    # no pedigree, no trio flag
    run_cohort_filter(tmp_path / "sites.vcf.gz", gt_paths, tmp_path / "c3.vcf.gz")
    out = {r.id: r for r in pysam.VariantFile(str(tmp_path / "c3.vcf.gz"))}
    assert all(list(r.filter.keys()) == ["PASS"] for r in out.values())
    assert "TR_OPP" not in out["G"].info


def test_weak_genotype():
    """No discordant pair and a copy ratio near 1: the de novo genotypes of
    dataset_paper children (DEL 1-10 kb PE 0 RD ~0.75-1.0, inherited PE ~10
    RD ~0.6); a deletion from 10 kb needs RD <= 0.62; below 1 kb never."""
    import numpy as np

    from svengers.cohort.metrics import weak_genotype

    pe = np.array([0, 0, 3, 0, 0])
    rd = np.array([0.9, 0.55, 0.9, np.nan, 0.7])
    assert weak_genotype("DEL", 5_000, pe, rd).tolist() == [True, False, False, False, False]
    assert weak_genotype("DEL", 50_000, pe, rd).tolist() == [True, False, False, False, True]
    assert not weak_genotype("DEL", 800, pe, rd).any()
    assert weak_genotype("DUP", 5_000, pe, np.array([1.1, 1.5, 1.1, np.nan, 1.3])).tolist() == [
        True, False, False, False, False]
    assert not weak_genotype("INV", 5_000, pe, rd).any()
    assert weak_genotype("DEL", 5_000, pe, rd, {"weak_del_rd": 0.95}).tolist() == [False] * 5
