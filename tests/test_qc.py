import numpy as np
import pytest

from vedette.qc import WAVE_MAX_K, contig_wave_score, sample_wave_score

BIN = 1000
N = 60_000  # 60 Mb contig at 1 kb bins


def _flat(rng, n=N):
    return rng.normal(1.0, 0.05, n)


def _wavy(rng, n=N, period_bins=800, amp=0.15):
    x = np.arange(n)
    return 1.0 + amp * np.sin(2 * np.pi * x / period_bins) + rng.normal(0, 0.05, n)


def test_flat_profile_scores_low():
    rng = np.random.default_rng(0)
    s = contig_wave_score(_flat(rng), np.zeros(N, bool), BIN)
    # white noise: energy fraction ~= band width / spectrum width ~ 1%
    assert s is not None and s < 0.10


def test_wavy_profile_scores_high():
    rng = np.random.default_rng(0)
    s = contig_wave_score(_wavy(rng), np.zeros(N, bool), BIN)  # 800 kb period
    assert s is not None and s > 0.5


def test_wave_score_orders_amplitudes():
    rng = np.random.default_rng(1)
    weak = contig_wave_score(_wavy(rng, amp=0.05), np.zeros(N, bool), BIN)
    strong = contig_wave_score(_wavy(rng, amp=0.3), np.zeros(N, bool), BIN)
    assert weak < strong


def test_masked_stretches_do_not_fabricate_waves():
    rng = np.random.default_rng(2)
    ratio = _flat(rng)
    mask = np.zeros(N, bool)
    mask[::7] = True          # scattered masked bins
    mask[10_000:12_000] = True  # one long masked block (interp would ramp here)
    s = contig_wave_score(ratio, mask, BIN)
    assert s is not None and s < 0.10


def test_short_contig_returns_none():
    rng = np.random.default_rng(3)
    assert contig_wave_score(_flat(rng, 3_000), np.zeros(3_000, bool), BIN) is None


def test_cn_step_is_not_a_wave():
    # a real 5 Mb heterozygous deletion is a step, not periodicity: the linear
    # detrend won't remove it entirely but it must not reach the wavy range
    rng = np.random.default_rng(4)
    ratio = _flat(rng)
    ratio[20_000:25_000] *= 0.5
    s = contig_wave_score(ratio, np.zeros(N, bool), BIN)
    wavy = contig_wave_score(_wavy(rng), np.zeros(N, bool), BIN)
    assert s < wavy


def test_sample_score_aggregates_and_gates_on_contig_count():
    rng = np.random.default_rng(5)
    ratio = {f"chr{i}": _wavy(rng) for i in range(1, 13)}
    mask = {k: np.zeros(N, bool) for k in ratio}
    score, per_contig = sample_wave_score(ratio, mask, BIN)
    assert score is not None and score > 0.5 and len(per_contig) == 12
    one = {"chr20": ratio["chr1"]}
    score_one, per_one = sample_wave_score(one, {"chr20": mask["chr1"]}, BIN)
    assert score_one is None and "chr20" in per_one


def test_band_definition_sane():
    assert WAVE_MAX_K == 50


def test_impact_labels_parse_as_scores(tmp_path):
    # every impact label must start with its score (the de novo tables sort
    # by that prefix; a bare label crashed the cohort report)
    import pysam

    from vedette.report import _GeneIndex, _impact_points

    gff = tmp_path / "genes.gff"
    lines = [
        'chr1\thavana\tgene\t10000\t20000\t.\t+\t.\tgene_id="ENSG1";gene_name="GENE1";gene_biotype="protein_coding";',
        'chr1\thavana\texon\t10000\t11000\t.\t+\t.\tgene_id="ENSG1";gene_name="GENE1";gene_biotype="protein_coding";',
    ]
    gff.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(gff), preset="gff", force=True)  # bgzips to .gff.gz
    idx = _GeneIndex(str(gff) + ".gz")
    for span in ((10_500, 10_600), (15_000, 16_000), (50_000, 60_000)):
        for svtype in ("DEL", "DUP"):
            score, label = idx.impact("chr1", *span, svtype)
            assert _impact_points(label) == score


def _genetrek_fixture(tmp_path):
    """A GFF with five coding genes and a GeneTrek table describing them."""
    import pysam

    genes = [  # name, start, end
        ("SHANK3", 100_000, 160_000), ("RECG", 200_000, 220_000),
        ("CANDG", 300_000, 320_000), ("CONSTR", 400_000, 420_000),
        ("PLAIN", 500_000, 520_000),
    ]
    lines = []
    for i, (name, a, b) in enumerate(genes):
        attr = f'gene_id="ENSG{i}";gene_name="{name}";gene_biotype="protein_coding";'
        lines.append(f"chr22\thavana\tgene\t{a}\t{b}\t.\t+\t.\t{attr}")
        lines.append(f"chr22\thavana\texon\t{a}\t{a + 1000}\t.\t+\t.\t{attr}")
        lines.append(f"chr22\thavana\texon\t{b - 1000}\t{b}\t.\t+\t.\t{attr}")
    gff = tmp_path / "genes.gff"
    gff.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(gff), preset="gff", force=True)
    cols = ["ensg", "hgnc_symbol", "hc_ndd_v7", "sysndd_definitive_v0_11_2",
            "candidate_ndd_v5", "inheritance", "shet_mean", "loeuf", "is_shet_constrained"]
    rows = [
        ["ENSG0", "SHANK3", "true", "true", "false", "Autosomal Dominant", "0.78", "0.37", "true"],
        ["ENSG1", "RECG", "true", "false", "false", "Autosomal Recessive", "0.01", "1.1", "false"],
        ["ENSG2", "CANDG", "false", "false", "true", "", "0.2", "0.5", "true"],
        ["ENSG3", "CONSTR", "false", "false", "false", "", "0.3", "0.2", "true"],
        ["ENSG4", "PLAIN", "false", "false", "false", "", "0.01", "1.5", "false"],
    ]
    tsv = tmp_path / "genetrek.tsv"
    tsv.write_text("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return str(gff) + ".gz", str(tsv)


def test_impact_score_follows_the_clingen_grid(tmp_path):
    """The old 0-4 score put half the coding genome on one level; the new one
    separates DEL from DUP and reads dosage and inheritance from GeneTrek."""
    from vedette.report import _GeneIndex

    gff, tsv = _genetrek_fixture(tmp_path)
    idx = _GeneIndex(gff, tsv)
    def pts(a, b, t="DEL"):
        return idx.impact("chr22", a, b, t)[0]
    assert pts(90_000, 170_000) == 1.00            # SHANK3 fully deleted
    assert pts(90_000, 120_000) == 0.90            # its first exon
    assert pts(120_000, 130_000) == 0.0            # intronic only
    assert pts(190_000, 230_000) == 0.15           # recessive NDD gene
    assert pts(290_000, 330_000) == 0.45           # candidate NDD, constrained
    assert pts(390_000, 430_000) == 0.15           # constrained, no NDD evidence
    assert pts(490_000, 530_000) == 0.0            # plain gene
    assert pts(600_000, 700_000) == -0.60          # no coding gene
    assert pts(100_500, 159_500, "DUP") == 0.45    # intragenic DUP of SHANK3
    assert pts(90_000, 170_000, "DUP") == 0.15     # whole-gene DUP: TS unknown
    label = idx.impact("chr22", 90_000, 170_000)[1]
    assert label.startswith("+1.00 — élevé") and "SHANK3" in label and "AD" in label


def _clingen_fixture(tmp_path):
    """ClinGen lists in the published layout: two comment lines, then a
    '#'-prefixed header."""
    reg = tmp_path / "regions.tsv"
    reg.write_text(
        "#ClinGen Region Curation Results\n#Genomic Locations are reported on GRCh38\n"
        "#ISCA ID\tISCA Region Name\tcytoBand\tGenomic Location\tHaploinsufficiency Score\t"
        "Triplosensitivity Score\n"
        "ISCA-1\t22q99 recurrent region (includes PLAIN)\t22q99\tchr22:1000001-2000000\t3\t3\n"
        "ISCA-2\t22q98 population region\t22q98\tchr22:3000001-3500000\t40\t0\n"
        "ISCA-3\tpending region\t22q97\ttbd\t3\t3\n"
    )
    genes = tmp_path / "genes.tsv"
    genes.write_text(
        "#ClinGen Gene Curation Results\n"
        "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\t"
        "Triplosensitivity Score\n"
        "PLAIN\t1\t22q\tchr22:500001-520000\t3\t3\n"
        "RECG\t2\t22q\tchr22:200001-220000\t30\tNot yet evaluated\n"
    )
    return str(genes), str(reg)


def test_impact_score_reads_the_clingen_dosage_map(tmp_path):
    """A recurrent region and ClinGen gene scores: the v5 review found a de
    novo 15q11.2q13 duplication (TS 3 region) scored +0.15."""
    from vedette.clingen import ClinGen
    from vedette.export import impact_fields
    from vedette.report import _GeneIndex

    gff, tsv = _genetrek_fixture(tmp_path)
    cg = ClinGen.load(*_clingen_fixture(tmp_path))
    assert cg.n_regions == 2 and cg.gene("RECG") == (30, None)
    idx = _GeneIndex(gff, tsv, cg)
    def imp(a, b, t="DEL"):
        return idx.impact("chr22", a, b, t)
    score, label = imp(1_020_000, 1_980_000, "DUP")   # 96 % of the TS 3 region
    assert score == 1.00 and "DUP complète de la région ClinGen 22q99 recurrent region" in label
    assert "(ISCA-1, TS 3)" in label and "includes" not in label
    assert imp(1_000_000, 1_300_000, "DUP")[0] == 0.45    # 30 %: partial, kept visible
    assert imp(1_000_000, 1_100_000, "DUP")[0] == -0.60   # 10 %, no gene: nothing
    assert imp(3_100_000, 3_200_000) == (-1.00, "-1.00 — bénin · dans la région ClinGen "
                                         "22q98 population region (ISCA-2, HI 40)")
    assert imp(3_100_000, 3_200_000, "DUP")[0] == -0.60   # benign for losses only
    assert impact_fields(imp(3_100_000, 3_200_000)[1])[:2] == ("-1.00", "benign")
    # ClinGen genes: PLAIN has no GeneTrek evidence, HI 3 and TS 3 in ClinGen
    score, label = imp(490_000, 530_000)
    assert score == 1.00 and "DEL complète PLAIN (ClinGen HI 3)" in label
    assert imp(490_000, 530_000, "DUP")[0] == 1.00
    assert imp(490_000, 505_000)[0] == 0.90                # exonic
    assert imp(190_000, 230_000)[0] == 0.15                # RECG: AR either way
    # without ClinGen, unchanged
    assert _GeneIndex(gff, tsv).impact("chr22", 490_000, 530_000)[0] == 0.0


def test_denovo_outlier_guard():
    from vedette.report import denovo_outliers

    def dn_rows(child, n, contig="chr1", step=1_000_000):
        return [{"child": child, "contig": contig, "svtype": "DEL",
                 "start": i * step, "end": i * step + 1_000} for i in range(n)]

    dn = {
        "F1": dn_rows("kidA", 3),
        "F2": dn_rows("kidB", 25) + dn_rows("kidC", 1),
        "F3": dn_rows("kidB", 5, contig="chr2"),   # counted across families
    }
    out = denovo_outliers(dn, threshold=20)
    assert out == {"kidB": 30}
    assert denovo_outliers({}, 20) == {}


def test_denovo_guard_counts_events_not_sites():
    """A de novo 22q13 deletion spans 50-80 cohort sites; counted as sites it
    made 77 PMS probands look like technical outliers."""
    from vedette.report import denovo_event_count, denovo_outliers

    # 60 nested/overlapping sites of one 5 Mb terminal deletion
    big = [{"child": "kid", "contig": "chr22", "svtype": "DEL",
            "start": 45_146_000 + i * 80_000, "end": 50_800_002} for i in range(60)]
    assert denovo_event_count(big) == 1
    assert denovo_outliers({"F": big}, threshold=20) == {}
    # the same intervals as DUPs are a second event, not the same one
    dups = [dict(r, svtype="DUP") for r in big[:3]]
    assert denovo_event_count(big + dups) == 2


def test_drop_karyotype_denovo():
    """An XXY boy's extra X reads as dozens of de novo chrX gains: they leave
    the de novo candidates, his autosomal ones stay; XYY drops chrY only."""
    from vedette.report import denovo_outliers, drop_karyotype_denovo

    def dn_rows(child, n, contig, step=1_000_000):
        return [{"child": child, "contig": contig, "svtype": "DUP",
                 "start": i * step, "end": i * step + 900_000} for i in range(n)]

    dn = {"F1": dn_rows("xxy", 40, "chrX") + dn_rows("xxy", 2, "chr7"),
          "F2": dn_rows("xyy", 3, "chrY") + dn_rows("xyy", 1, "chrX"),
          "F3": dn_rows("xy", 2, "chrX")}
    karyo = {"xxy": {"karyotype": "XXY"}, "xyy": {"karyotype": "XYY"},
             "xy": {"karyotype": "XY"}}
    assert denovo_outliers(dn, threshold=20) == {"xxy": 42}
    removed = drop_karyotype_denovo(dn, karyo)
    assert removed == {"xxy": ("XXY", 40), "xyy": ("XYY", 3)}
    assert [r["contig"] for r in dn["F1"]] == ["chr7", "chr7"]
    assert [r["contig"] for r in dn["F2"]] == ["chrX"]
    assert len(dn["F3"]) == 2
    assert denovo_outliers(dn, threshold=20) == {}


def test_family_events_page(tmp_path):
    from vedette.report import _write_events_page

    members = ["C000EZH", "C000EZI", "C000EZJ"]
    ev = {
        "contig": "chr1", "pos": 1000, "end": 2000, "svtype": "DEL", "size": 1000,
        "filter": "PASS", "rescue": "", "genes": "GENE1", "impact": "3 — exonique, NDD-candidat",
        "inherit": "de novo", "carriers": ["C000EZJ"], "sqmax": "12.3", "cfrac": "0.333",
        "refcf": "0.004", "lcr": "0.00", "segdup": "1.00", "own": "PASS", "support": "fort",
        "gt": [("0/0", 99, "2"), ("0/0", 99, "2"), ("0/1", 45, "1")],
    }
    ev2 = dict(ev, contig="chr2", svtype="DUP", filter="LowSQ", carriers=[],
               rescue="RECURRENT", inherit="maternel")
    _write_events_page(tmp_path / "events.html", "Famille test", members, [ev, ev2])

    page = (tmp_path / "events.html").read_text()
    assert "1 DEL + 1 DUP" in page
    # the chr/type selectors and the rows' filter attributes drive the JS
    assert '<select id="chr">' in page and '<select id="type">' in page
    assert '<option value="chr1">chr1 (1)</option>' in page
    assert 'data-chr="chr2"' in page and 'data-type="DUP" data-pass="0"' in page
    assert "GENE1" in page and "NDD-candidat" in page and "RECURRENT" in page and "nopass" in page
    # inheritance column + clickable-header sorting
    assert "<th>héritage</th>" in page and "de novo" in page and "maternel" in page
    assert "th.onclick" in page and "localeCompare" in page
    # header outside tbody: the filter JS must never hide or sort it
    assert "<thead>" in page and "</thead><tbody>" in page
    # manual validation + IGV: per-row controls (vrai/faux/incertain), export/import in the bar
    assert 'class="v-true"' in page and 'class="v-false"' in page and 'class="v-unsure"' in page
    assert 'class="igv"' in page and 'option value="unsure"' in page
    assert 'id="export"' in page and 'id="import"' in page and 'data-site="chr1:1000-2000"' in page
    # cohort-level context columns (merge SQ, run carrier fraction, reference-cohort frequency)
    assert "<th>SQmax</th>" in page and "<th>% port.</th>" in page and "<th>fréq. réf.</th>" in page
    assert "<th>LCR</th>" in page and "<th>segdup</th>" in page
    assert 'data-key="support">support' in page and 'data-support="fort"' in page
    assert "<td>12.3</td><td>0.333</td><td>0.004</td><td>0.00</td><td>1.00</td>" in page
    assert "GT<span" in page
    # the TSV export resolves columns by header text, not by cell index
    assert 'COL[th.textContent.trim()]' in page and '["lcr", "LCR"]' in page


def test_carrier_excess_guard():
    """The guard must key on carrier LOAD, not on the low-GQ fraction: on
    dataset_paper the fraction's cohort median is 0.75 and the three real
    offenders sit at 0.61-0.69 (percentiles 2-6), so a fraction threshold
    flagged 2670/2670 samples."""
    import numpy as np

    from vedette.report import carrier_excess

    n_carrier = np.array([300, 310, 305, 1830, 295])   # one sample at 6x
    exc = carrier_excess(n_carrier)
    assert exc[0] == pytest.approx(300 / 305)
    assert exc[3] == pytest.approx(1830 / 305)
    assert (exc >= 5.0).sum() == 1        # production threshold
    assert (exc >= 1.3).sum() == 1

    # degenerate cohort (no carriers at all) -> all NaN, never a flag
    assert np.isnan(carrier_excess(np.array([0, 0]))).all()


def test_low_gq_carrier_fraction():
    import numpy as np

    from vedette.report import pass_carrier_stats

    # 4 sites x 3 samples; site 2 is non-PASS (must be excluded from the calc)
    nalt = np.array([
        [1, 0, -1],   # site0 PASS: S0 het carrier, S1 ref, S2 missing
        [2, 1, 0],    # site1 PASS: all carriers except S2
        [2, 2, 2],    # site2 NOT PASS: would look terrible if counted
        [0, 1, 1],    # site3 PASS: S1/S2 carriers
    ], dtype=np.int8)
    gq = np.array([
        [5, 30, 0],
        [10, 40, 30],
        [1, 1, 1],
        [30, 15, 45],
    ], dtype=np.int16)
    pass_mask = np.array([True, True, False, True])

    n_carrier, frac = pass_carrier_stats(nalt, gq, pass_mask)
    assert list(n_carrier) == [2, 2, 1]
    # S0 PASS carriers: site0(GQ5<20), site1(GQ10<20) -> 2/2 = 1.0
    assert frac[0] == 1.0
    # S1 PASS carriers: site1(GQ40 ok), site3(GQ15<20) -> 1/2 = 0.5
    assert frac[1] == 0.5
    # S2 PASS carriers: site3(GQ45 ok) only (site0/site1 not a carrier) -> 0/1 = 0.0
    assert frac[2] == 0.0

    # a sample with zero PASS carriers is NaN, not zero (nothing to judge)
    nalt_none = np.array([[0], [0]], dtype=np.int8)
    gq_none = np.array([[10], [10]], dtype=np.int16)
    _, frac_none = pass_carrier_stats(nalt_none, gq_none, np.array([True, True]))
    assert np.isnan(frac_none[0])


def test_collapse_overlapping_events_per_carrier_set():
    """One biological deletion shows up as a stack of nested cohort sites when
    patients have different breakpoints (PMS 22q13). Merge them — but only
    when the carrier set is identical, so an inherited event never absorbs a
    de novo one it overlaps."""
    from vedette.report import collapse_events

    def ev(pos, end, carriers, sid, **kw):
        return dict(id=sid, contig="chr22", pos=pos, end=end, svtype="DEL",
                    size=end - pos, filter="PASS", rescue="", genes="SHANK3",
                    impact="4 — exonique", inherit=kw.get("inherit", ""),
                    sqmax=kw.get("sqmax", ""), cfrac="", refcf="", lcr="", segdup="",
                    carriers=carriers, gt=[])

    child = ["SMP-08"]
    trio = ["SMP-08", "SMP-09", "SMP-10"]
    out = collapse_events([
        ev(50_676_403, 50_757_929, child, "s1", sqmax="8.0"),
        ev(50_686_127, 50_719_326, child, "s2", sqmax="10.9"),
        ev(50_722_582, 50_736_717, child, "s3", sqmax="17.2"),
        ev(50_692_669, 50_692_985, trio, "s4"),       # biparental: stays apart
        ev(51_000_000, 51_000_500, child, "s5"),      # disjoint: stays apart
    ])
    merged = {(e["pos"], e["end"]): e for e in out}
    assert len(out) == 3
    big = merged[(50_676_403, 50_757_929)]
    assert big["n_sites"] == 3 and set(big["_ids"]) == {"s1", "s2", "s3"}
    assert big["size"] == 50_757_929 - 50_676_403
    assert big["sqmax"] == "17.2"                     # strongest evidence kept
    assert merged[(50_692_669, 50_692_985)]["n_sites"] == 1   # different carriers
    assert merged[(51_000_000, 51_000_500)]["n_sites"] == 1   # no overlap


def test_own_call_status():
    """A cohort site is a breakpoint cluster; re-genotyping then genotypes
    every sample, backed by an own call or not. On PMS 80% of a family's PASS
    events have no own call and blind IGV review scored those 0/33."""
    from vedette.report import own_call_status

    idx = {
        ("chr22", "DEL"): [(50_000, 60_000, False), (80_000, 90_000, True)],
        ("chr22", "DUP"): [(50_000, 60_000, True)],
    }
    assert own_call_status(idx, "chr22", 80_100, 89_900, "DEL") == "PASS"
    assert own_call_status(idx, "chr22", 50_100, 59_900, "DEL") == "filtré"
    assert own_call_status(idx, "chr22", 50_100, 59_900, "DUP") == "PASS"  # type matters
    assert own_call_status(idx, "chr22", 200_000, 210_000, "DEL") == ""
    # reciprocal overlap below 0.5 does not count as backing
    assert own_call_status(idx, "chr22", 59_000, 120_000, "DEL") == ""
    assert own_call_status(idx, "chr1", 50_100, 59_900, "DEL") == ""


def test_own_span_takes_the_carriers_breakpoints():
    """A cohort site carries its cluster representative's coordinates. PMS
    pools 22q13 deletions of every size, so FAM-10's own call
    (45,146,001; IGV 45,145,667) was shown at its site's 44,116,001."""
    from vedette.report import own_span

    idx = {("chr22", "DEL"): [
        (44_200_000, 48_900_000, True, "mum"),    # not a carrier of the row
        (44_617_439, 44_617_645, True, "kid"),    # small polymorphism: not chained
        (45_146_001, 48_986_001, True, "kid"),
    ]}
    assert own_span(idx, "chr22", 44_116_001, 48_986_002, "DEL", {"kid"}) == (45_146_001, 48_986_001)
    assert own_span(idx, "chr22", 44_116_001, 48_986_002, "DEL", {"dad"}) is None
    # filtered own calls never move a row
    idx2 = {("chr22", "DEL"): [(1_000, 9_000, False, "kid")]}
    assert own_span(idx2, "chr22", 1_000, 10_000, "DEL", {"kid"}) is None
    # own evidence under half of the site: the site's bounds stay
    idx3 = {("chr22", "DEL"): [(1_000, 4_000, True, "kid")]}
    assert own_span(idx3, "chr22", 1_000, 10_000, "DEL", {"kid"}) is None


def test_own_span_chains_a_fragmented_own_call():
    """FAM-09 (IGV 43,784,666): the proband's deletion comes back as
    nine own PASS calls from 43,785,001, each overlapping the 43,570,001 site
    by less than half — chained, they match it."""
    from vedette.report import own_span

    idx = {("chr22", "DEL"): [
        (43_153_340, 43_153_667, True, "kid"),    # unrelated, before the site
        (43_785_001, 44_113_001, True, "kid"),
        (44_116_001, 44_739_001, True, "kid"),
        (44_742_001, 48_988_001, True, "kid"),
    ]}
    assert own_span(idx, "chr22", 43_570_001, 46_859_002, "DEL", {"kid"}) == (43_785_001, 48_988_001)


def test_own_span_grows_to_the_carriers_breakpoint_but_keeps_scale():
    """FAM-11's duplication was merged into a site starting 141 kb
    inside it: the row must reach the carrier's own breakpoint, or the event
    splits in two. But a small site inside a large own call — here a
    polymorphism the mother carries too — keeps its own bounds."""
    from vedette.report import own_span

    idx = {("chr22", "DUP"): [(43_835_001, 46_859_001, True, "kid")]}
    assert own_span(idx, "chr22", 43_975_501, 46_732_502, "DUP", {"kid"}) == (43_835_001, 46_859_001)
    assert own_span(idx, "chr22", 45_000_000, 45_000_700, "DUP", {"kid", "mum"}) is None


def test_own_span_links_through_filtered_pieces_but_not_past_the_edge():
    """FAM-08: the middle 138 kb of the deletion is a LowSQ own call;
    it links the two PASS pieces but never sets an end. And a subtelomeric
    PASS call sharing 1 bp with the site must not push the row past chr22."""
    from vedette.report import own_span

    idx = {("chr22", "DEL"): [
        (46_069_001, 47_611_001, True, "kid"),
        (47_612_001, 47_750_001, False, "kid"),
        (47_753_001, 48_986_001, True, "kid"),
    ]}
    assert own_span(idx, "chr22", 46_069_001, 48_986_002, "DEL", {"kid"}) == (46_069_001, 48_986_001)
    # a chain of filtered calls only: nothing to anchor the ends
    only_filtered = {("chr22", "DEL"): [(1_000, 9_000, False, "kid")]}
    assert own_span(only_filtered, "chr22", 1_000, 9_000, "DEL", {"kid"}) is None
    tel = {("chr22", "DEL"): [(50_346_001, 50_800_001, True, "kid"),
                              (50_800_001, 50_819_001, True, "kid")]}
    assert own_span(tel, "chr22", 50_346_001, 50_800_002, "DEL", {"kid"}) == (50_346_001, 50_800_001)


def test_own_span_ignores_own_calls_of_failed_cohort_sites():
    """FAM-14: a 23 kb own PASS call in the subtelomere belongs to the
    site every PMS sample carries (non-PASS at cohort level); it must not set
    the row's end past the end of chr22."""
    from vedette.report import _cohort_site_pass, _site_index, own_span

    sites = [
        {"contig": "chr22", "svtype": "DEL", "start": 50_346_000, "end": 50_800_002, "filter": "PASS"},
        {"contig": "chr22", "svtype": "DEL", "start": 50_800_000, "end": 50_819_002,
         "filter": "CohortArtifact;HighCarrierFrac"},
    ]
    sidx = _site_index(sites)
    assert _cohort_site_pass(sidx, "chr22", 50_346_001, 50_796_001, "DEL") is True
    assert _cohort_site_pass(sidx, "chr22", 50_796_001, 50_819_001, "DEL") is False
    assert _cohort_site_pass(sidx, "chr22", 1_000, 2_000, "DEL") is None
    idx = {("chr22", "DEL"): [
        (50_346_001, 50_796_001, True, "kid", True),
        (50_796_001, 50_819_001, True, "kid", False),
    ]}
    assert own_span(idx, "chr22", 50_346_001, 50_800_002, "DEL", {"kid"}) == (50_346_001, 50_796_001)


def test_own_span_nested_calls_do_not_set_ends():
    """FAM-06: small PASS calls nested in a large LowSQ duplication
    must not become the row's breakpoints (the 5.6 Mb duplication split)."""
    from vedette.report import own_span

    idx = {("chr8", "DUP"): [
        (140_623_001, 141_335_001, False, "kid", None),
        (140_623_177, 140_623_293, True, "kid", True),
        (141_091_418, 141_091_983, True, "kid", True),
        (141_308_477, 141_308_784, True, "kid", True),
    ]}
    assert own_span(idx, "chr8", 140_623_001, 141_335_002, "DUP", {"kid"}) is None


def test_own_call_index_clamps_ends_to_the_contig(tmp_path):
    """Per-sample VCFs written before vcfout clamped END overhang the contig
    by up to one 1 kb bin (dataset_paper: chrX 156,041,002 for 156,040,895);
    rows take own-call ends directly, so the index must bound them."""
    import pysam

    from vedette.report import own_call_index

    d = tmp_path / "samples" / "kid"
    d.mkdir(parents=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chrX", length=156_040_895)
    hdr.add_line('##ALT=<ID=DEL,Description="Deletion">')
    hdr.info.add("SVTYPE", "1", "String", "type")
    hdr.info.add("END", "1", "Integer", "end")
    hdr.add_sample("kid")
    with pysam.VariantFile(str(d / "kid.filtered.vcf.gz"), "wz", header=hdr) as vf:
        rec = vf.new_record(contig="chrX", start=156_010_000, stop=156_041_001,
                            alleles=("N", "<DEL>"))
        rec.info["SVTYPE"] = "DEL"
        rec.filter.add("PASS")
        vf.write(rec)
    idx = own_call_index(tmp_path, ["kid"])
    ((pos, end, is_pass, who),) = idx[("chrX", "DEL")]
    assert end == 156_040_895 and is_pass and who == "kid"


def test_read_cohort_vcf_clamps_site_ends_to_the_contig(tmp_path):
    """Cohort sites merged before merge clamped END overhang the contig
    (dataset_paper: 24 sites, chr12 133,276,002 for 133,275,309); a family
    row without an own span keeps the site's ends, so the reader bounds them."""
    import pysam

    from vedette.report import _read_cohort_vcf

    chr12 = 133_275_309
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chr12", length=chr12)
    hdr.add_line('##ALT=<ID=DEL,Description="Deletion">')
    hdr.info.add("SVTYPE", "1", "String", "type")
    hdr.info.add("END", "1", "Integer", "end")
    for key, typ in (("AF", "Float"), ("NS", "Integer"), ("NCARRIER", "Integer"),
                     ("HWE_P", "Float"), ("RESCUE", "String")):
        hdr.info.add(key, "1", typ, key)
    hdr.formats.add("GT", "1", "String", "Genotype")
    hdr.formats.add("GQ", "1", "Integer", "Genotype quality")
    hdr.add_sample("kid")
    path = tmp_path / "cohort.vcf.gz"
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for start, stop in ((1_000_000, 1_020_001), (133_250_000, 133_276_002)):
            rec = vf.new_record(contig="chr12", start=start, stop=stop,
                                alleles=("N", "<DEL>"))
            rec.info["SVTYPE"] = "DEL"
            rec.samples["kid"]["GT"] = (0, 1)
            rec.samples["kid"]["GQ"] = 99
            vf.write(rec)
    _, sites, *_ = _read_cohort_vcf(path)
    by_start = {st["start"]: st for st in sites}
    assert by_start[133_250_000]["end"] == chr12
    assert by_start[133_250_000]["size"] == chr12 - 133_250_000
    assert by_start[1_000_000]["end"] == 1_020_000  # inside: untouched


def test_sample_qc_flags_and_family_row_alert(tmp_path):
    """dataset_paper: family FAM-03 held a dozen of the 57 rows at
    impact >= +1.00, its children carrying 478 and 94 de novo candidates. The
    row must carry its carriers' QC alerts, filterable on the page."""
    import numpy as np

    from vedette.report import _FamilySpill, _write_all_family_events, sample_qc_flags

    flags = sample_qc_flags({"kid": 12.3}, {"kid": 478})
    assert flags == {"kid": ["charge 12.3×", "de novo 478"]}

    samples = ["kid", "dad", "mum"]
    sites = [
        {"id": "s1", "contig": "chr1", "start": 999, "end": 2000, "svtype": "DEL", "size": 1001,
         "filter": "PASS", "rescue": "", "idx": 0, "sqmax": None, "carrier_frac": None,
         "ref_carrier_frac": None},
        {"id": "s2", "contig": "chr1", "start": 9999, "end": 12000, "svtype": "DEL", "size": 2001,
         "filter": "PASS", "rescue": "", "idx": 1, "sqmax": None, "carrier_frac": None,
         "ref_carrier_frac": None},
    ]
    nalt = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.int8)     # s1 kid, s2 dad
    gq = np.array([[60, 99, 99], [99, 60, 99]], dtype=np.int16)
    cn = np.array([[1, 2, 2], [2, 1, 2]], dtype=np.int8)
    wave = {"kid": 0.2, "dad": 0.05, "mum": None}
    n = _write_all_family_events(tmp_path, {"F1": samples}, samples, sites, nalt, gq, cn,
                                 lambda *a: "", lambda *a: (0, ""), _FamilySpill(tmp_path / "inh"),
                                 sample_qc=flags, wave_scores=wave)
    assert n == 1
    page = (tmp_path / "families" / "F1" / "events.html").read_text()
    assert page.count('data-qc="alerte"') == 1 and page.count('data-qc="ok"') == 1
    assert "⚠ kid wave 0.200, charge 12.3×, de novo 478" in page
    assert '<select id="qc">' in page


def test_row_wave_score_leaves_out_the_rows_contig():
    """PMS: a >= 1 Mb chr22 deletion alone pushed 26 of 48 carriers over the
    wave threshold; without chr22 in the mean, 5 remain."""
    from vedette.qc import WAVE_WARN
    from vedette.report import row_wave_score

    flat = {f"chr{i}": 0.03 for i in range(1, 22)}
    proband = flat | {"chr22": 0.49, "chr17": 0.58, "chr16": 0.37}   # deletion + GC-rich
    assert row_wave_score(proband, 0.116, "chr22") < row_wave_score(proband, 0.116, "chr1")
    wavy = {f"chr{i}": 0.30 for i in range(1, 23)}
    assert row_wave_score(wavy, 0.30, "chr22") >= WAVE_WARN       # a wavy genome stays wavy
    assert row_wave_score({}, 0.2, "chr22") == 0.2               # old stores: sample score
    assert row_wave_score(None, None, "chr22") is None


def test_family_page_uses_the_carriers_own_breakpoints(tmp_path):
    """End to end: the row of a re-genotyped site backed by the child's own
    call is drawn at the own call's coordinates, with its genes recomputed."""
    import numpy as np
    import pysam

    from vedette.report import _FamilySpill, _write_all_family_events

    samples = ["kid", "dad", "mum"]
    d = tmp_path / "samples" / "kid"
    d.mkdir(parents=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chr22", length=50_818_468)
    hdr.add_line('##ALT=<ID=DEL,Description="Deletion">')
    hdr.info.add("SVTYPE", "1", "String", "type")
    hdr.info.add("END", "1", "Integer", "end")
    hdr.add_sample("kid")
    with pysam.VariantFile(str(d / "kid.filtered.vcf.gz"), "wz", header=hdr) as vf:
        rec = vf.new_record(contig="chr22", start=45_146_000, stop=48_986_001,
                            alleles=("N", "<DEL>"))
        rec.info["SVTYPE"] = "DEL"
        rec.filter.add("PASS")
        vf.write(rec)
    sites = [{"id": "s1", "contig": "chr22", "start": 44_116_000, "end": 48_986_002,
              "svtype": "DEL", "size": 4_870_002, "filter": "PASS", "rescue": "",
              "idx": 0, "sqmax": None, "carrier_frac": None, "ref_carrier_frac": None}]
    nalt = np.array([[1, 0, 0]], dtype=np.int8)
    gq = np.array([[99, 99, 99]], dtype=np.int16)
    cn = np.array([[1, 2, 2]], dtype=np.int8)
    seen = []

    def genes(contig, start0, end):
        seen.append((start0, end))
        return "G"

    n = _write_all_family_events(tmp_path, {"F1": samples}, samples, sites, nalt, gq, cn,
                                 genes, lambda *a: (0, ""), _FamilySpill(tmp_path / "inh"))
    assert n == 1
    page = (tmp_path / "families" / "F1" / "events.html").read_text()
    assert "<td>45,146,001</td><td>48,986,000</td>" in page
    assert "44,116,001" not in page
    assert (45_146_000, 48_986_000) in seen      # genes of the own span


def test_carrier_requires_sufficient_gq():
    """Re-genotyping emits GT>0 with a low GQ far too readily: on PMS the
    apparent de novo rate — a genotype-error proxy — is 16.7% for DEL
    carriers below GQ 20 against 0.8% above (13.9% vs 1.9% for DUP), so a
    low-GQ alternative genotype must not make the sample a carrier."""
    import numpy as np

    from vedette.report import _FamilySpill, _write_all_family_events

    samples = ["kid", "dad", "mum"]
    sites = [{"id": "s1", "contig": "chr1", "start": 999, "end": 2000, "svtype": "DEL",
              "size": 1001, "filter": "PASS", "rescue": "", "idx": 0,
              "sqmax": None, "carrier_frac": None, "ref_carrier_frac": None}]
    #                kid carrier GQ 45, dad carrier GQ 3 (noise), mum ref
    nalt = np.array([[1, 1, 0]], dtype=np.int8)
    gq = np.array([[45, 3, 99]], dtype=np.int16)
    cn = np.array([[1, 1, 2]], dtype=np.int8)

    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp)
        spill = _FamilySpill(run / "inh")
        n = _write_all_family_events(
            run, {"F1": samples}, samples, sites, nalt, gq, cn,
            lambda *a: "", lambda *a: (0, ""), spill,
        )
        assert n == 1
        page = (run / "families" / "F1" / "events.html").read_text()
    # dad's low-GQ genotype is shown but he is not listed as a carrier
    assert "0/1<span" in page          # genotypes still displayed
    assert ">kid<" in page or "kid" in page
    import re
    carriers = re.search(r'<td class="c">([^<]*)</td><td><a class="igv"', page)
    assert carriers and "dad" not in carriers.group(1)


def test_large_deletion_keeps_pass_when_chromosome_depth_agrees(tmp_path):
    """Deletions past ~1 Mb are outside the model's training range (131 of
    132555 labelled DEL examples, all of them total-dropout regions), and once
    the HMM splits a terminal deletion into adjacent segments each segment's
    flanks fall inside its neighbour, so dhffc ~ 1 and the model sees nothing.
    fc_chrom still does, so it keeps the call."""

    from vedette.config import load_config

    cfg = load_config(None)
    f = cfg.filter
    # the real case: chr22:42,559,001-48,986,001 in SMP-07, SQ 3.73, CN=1
    big, small = 6_427_000, 120_000   # 120 kb: below the 200 kb rescue floor
    def override(size, fc_chrom, maskfrac):
        return (size >= f.large_del_rescue_size
                and fc_chrom < f.large_del_rescue_fc_chrom
                and maskfrac < f.large_del_rescue_max_mask)

    assert override(big, 0.528, 0.0)          # the rejected PMS deletion
    assert not override(small, 0.528, 0.0)    # too small: model is in range
    assert override(229_000, 0.528, 0.0)      # the 229 kb first segment, rescued at 200 kb
    assert not override(big, 0.98, 0.0)       # depth says no copy loss
    assert not override(big, 0.528, 0.8)      # mostly masked: untrustworthy


def test_collapse_bridges_small_gaps_but_not_distant_events():
    """A fragmented terminal deletion comes back as several intervals with
    CN=0 islands or filtered segments between them (FAM-01's 8.9 Mb
    22q13 deletion arrived as four pieces split by 3, 15 and 55 kb). The
    tolerated gap scales with the larger neighbour and is capped, so small
    events stay apart."""
    from vedette.report import GAP_ABS, GAP_FRAC, GAP_MAX, collapse_events

    def ev(pos, end, sid):
        return dict(id=sid, contig="chr22", pos=pos, end=end, svtype="DEL",
                    size=end - pos, filter="PASS", rescue="", genes="", impact="",
                    inherit="", sqmax="", cfrac="", refcf="", lcr="", segdup="",
                    carriers=["SMP-07"], gt=[])

    # 55 kb gap beside a 6.43 Mb segment: 2% of it is 129 kb, so bridged
    out = collapse_events([ev(42_138_001, 42_504_002, "b"), ev(42_559_001, 48_988_002, "c")])
    assert len(out) == 1 and out[0]["n_sites"] == 2
    assert (out[0]["pos"], out[0]["end"]) == (42_138_001, 48_988_002)

    # same 55 kb gap between two 366 kb segments: 2% is 7 kb, so kept apart
    apart = collapse_events([ev(42_138_001, 42_504_002, "b"), ev(42_559_001, 42_925_002, "d")])
    assert len(apart) == 2

    # the floor bridges sub-2 kb gaps whatever the size
    tiny = collapse_events([ev(1_000_000, 1_005_000, "x"), ev(1_006_500, 1_010_000, "y")])
    assert len(tiny) == 1

    # a small event does not bridge a gap longer than itself: FAM-10's
    # 0.6 kb polymorphism 21 kb before its 4.9 Mb deletion stayed glued to it
    lone = collapse_events([ev(44_094_543, 44_095_118, "p"), ev(44_116_001, 48_986_002, "q")])
    assert len(lone) == 2

    # and the cap bounds the worst case: 2% of 50 Mb would be 1 Mb
    assert min(max(GAP_ABS, int(GAP_FRAC * 50_000_000)), GAP_MAX) == GAP_MAX


def test_merged_row_shows_majority_inheritance():
    """A terminal deletion picks up a few sites a parent also carries;
    concatenating every label buried the fact that the event is de novo."""
    from vedette.report import collapse_events

    def ev(pos, end, sid, inherit):
        return dict(id=sid, contig="chr22", pos=pos, end=end, svtype="DEL",
                    size=end - pos, filter="PASS", rescue="", genes="", impact="",
                    inherit=inherit, sqmax="", cfrac="", refcf="", lcr="", segdup="",
                    carriers=["SMP-07"], gt=[])

    out = collapse_events([
        ev(1_000_000, 2_000_000, "a", "de novo"),
        ev(2_000_500, 3_000_000, "b", "de novo"),
        ev(3_000_500, 4_000_000, "c", "paternel"),
    ])
    assert len(out) == 1
    assert out[0]["inherit"] == "de novo (2/3)"

    # unanimous: no fraction, just the label
    clean = collapse_events([
        ev(1_000_000, 2_000_000, "a", "maternel"),
        ev(2_000_500, 3_000_000, "b", "maternel"),
    ])
    assert clean[0]["inherit"] == "maternel"


def test_denovo_cutoff_is_cohort_relative():
    """An absolute 20 flagged 664 of 665 dataset_paper children; 5x the
    cohort median flags the tail only, with the floor for family runs."""
    from vedette.report import denovo_cutoff

    cohort = {f"k{i}": 100 for i in range(76)} | {"noisy": 576}
    assert denovo_cutoff(cohort, 20) == 500.0
    assert {c for c, n in cohort.items() if n > denovo_cutoff(cohort, 20)} == {"noisy"}
    assert denovo_cutoff({"kid": 3, "kid2": 40}, 20) == 20.0     # too few children


def test_big_cnv_alert_on_family_rows(tmp_path):
    """dataset_paper SMP-02: 13 non-transmitted deletions >= 1 Mb, missed by
    the wave, carrier-load and de novo guards; >= 3 large own-PASS events
    flags the carrier on every row it carries."""
    import numpy as np

    from vedette.report import BIG_CNV_WARN, _FamilySpill, _write_all_family_events, big_cnv_counts

    def ev(pos, size, own="PASS", flt="PASS", car=("dad",)):
        return dict(contig="chr10", pos=pos, end=pos + size, filter=flt, own=own,
                    carriers=list(car))
    n = big_cnv_counts([ev(1, 2_000_000), ev(5_000_000, 1_500_000), ev(9_000_000, 900_000),
                        ev(12_000_000, 3_000_000, own="filtré"), ev(20_000_000, 2_000_000, flt="LowSQ")])
    assert n["dad"] == 2                      # 900 kb, filtered own call, non-PASS: not counted

    samples = ["kid", "dad", "mum"]
    sites = []
    for i in range(BIG_CNV_WARN):
        start = 10_000_000 * (i + 1)
        sites.append({"id": f"s{i}", "contig": "chr10", "start": start, "end": start + 1_200_000,
                      "svtype": "DEL", "size": 1_200_000, "filter": "PASS", "rescue": "",
                      "idx": i, "sqmax": None, "carrier_frac": None, "ref_carrier_frac": None})
    nalt = np.array([[0, 1, 0]] * BIG_CNV_WARN, dtype=np.int8)
    gq = np.array([[99, 99, 99]] * BIG_CNV_WARN, dtype=np.int16)
    cn = np.array([[2, 1, 2]] * BIG_CNV_WARN, dtype=np.int8)
    # dad's own PASS calls behind each site
    import pysam
    d = tmp_path / "samples" / "dad"
    d.mkdir(parents=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chr10", length=133_797_422)
    hdr.add_line('##ALT=<ID=DEL,Description="Deletion">')
    hdr.info.add("SVTYPE", "1", "String", "t")
    hdr.info.add("END", "1", "Integer", "e")
    hdr.add_sample("dad")
    with pysam.VariantFile(str(d / "dad.filtered.vcf.gz"), "wz", header=hdr) as vf:
        for st in sites:
            r = vf.new_record(contig="chr10", start=st["start"], stop=st["end"], alleles=("N", "<DEL>"))
            r.info["SVTYPE"] = "DEL"
            r.filter.add("PASS")
            vf.write(r)
    _write_all_family_events(tmp_path, {"F1": samples}, samples, sites, nalt, gq, cn,
                             lambda *a: "", lambda *a: (0, ""), _FamilySpill(tmp_path / "inh"))
    page = (tmp_path / "families" / "F1" / "events.html").read_text()
    assert page.count('data-qc="alerte"') == BIG_CNV_WARN
    assert f"dad grands CNV {BIG_CNV_WARN}" in page



def _sub_vcf(run_dir, sample, calls):
    """Per-sample filtered VCF with (start, stop, svtype, filters, cn, gq, sq)."""
    import pysam

    d = run_dir / "samples" / sample
    d.mkdir(parents=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chr22", length=50_818_468)
    for t in ("DEL", "DUP"):
        hdr.add_line(f'##ALT=<ID={t},Description="{t}">')
    hdr.add_line('##FILTER=<ID=LowSQ,Description="low SQ">')
    hdr.add_line('##FILTER=<ID=LowMap,Description="low mappability">')
    hdr.info.add("SVTYPE", "1", "String", "type")
    hdr.info.add("END", "1", "Integer", "end")
    hdr.formats.add("GT", "1", "String", "Genotype")
    hdr.formats.add("GQ", "1", "Integer", "GQ")
    hdr.formats.add("CN", "1", "Integer", "CN")
    hdr.formats.add("SQ", "1", "Float", "SQ")
    hdr.add_sample(sample)
    with pysam.VariantFile(str(d / f"{sample}.filtered.vcf.gz"), "wz", header=hdr) as vf:
        for start, stop, svtype, filters, cn, gq, sq in calls:
            rec = vf.new_record(contig="chr22", start=start, stop=stop,
                                alleles=("N", f"<{svtype}>"))
            rec.info["SVTYPE"] = svtype
            for f in filters:
                rec.filter.add(f)
            s = rec.samples[sample]
            s["GT"], s["GQ"], s["CN"], s["SQ"] = "0/1", gq, cn, sq
            vf.write(rec)


def test_subthreshold_ndd_call_is_a_row_of_its_family_page(tmp_path):
    """A PMS proband's 9 kb deletion of the 3' end of SHANK3 (CN 1, GQ 99,
    confirmed in IGV) failed the score alone (SQ 0.38): it must reach the
    family page as a sub-threshold row; the same call in an unrelated genome
    (an artefact), a call with another filter, a CN-2 call and a call in no
    NDD gene must not."""
    from vedette.report import _FamilySpill, _write_all_family_events, subthreshold_candidates

    shank3 = (50_731_000, 50_740_001, "DEL", ["LowSQ"], 1, 99, 0.38)
    _sub_vcf(tmp_path, "kid", [
        shank3,
        (50_600_000, 50_610_001, "DEL", ["LowSQ"], 1, 99, 0.5),       # recurrent
        (50_650_000, 50_660_001, "DEL", ["LowSQ", "LowMap"], 1, 99, 0.5),
        (50_700_000, 50_705_001, "DEL", ["LowSQ"], 2, 99, 0.5),       # CN 2
        (40_000_000, 40_010_001, "DEL", ["LowSQ"], 1, 99, 0.5),       # no NDD gene
    ])
    _sub_vcf(tmp_path, "dad", [])
    _sub_vcf(tmp_path, "mum", [])
    _sub_vcf(tmp_path, "stranger", [(50_600_500, 50_610_001, "DEL", ["LowSQ"], 1, 95, 0.2)])

    samples = ["kid", "dad", "mum", "stranger"]
    cands = subthreshold_candidates(tmp_path, samples)
    assert {(c["start"], c["cn"]) for c in cands.calls("kid")} == {
        (50_731_000, 1), (50_600_000, 1), (40_000_000, 1)}

    def impact(contig, start, end, svtype="DEL"):
        if start < 50_733_298 and end > 50_674_415:
            return 0.9, "+0.90 — élevé · DEL exonique SHANK3"
        if start < 50_620_000 and end > 50_590_000:
            return 0.9, "+0.90 — élevé · DEL exonique XYZ"
        return 0.0, "+0.00 — faible"

    nalt = np.zeros((0, 4), dtype=np.int8)
    assert cands.sharers("kid", "chr22", 50_600_000, 50_610_001, "DEL") == {"stranger"}
    assert cands.sharers("kid", "chr22", 50_731_000, 50_740_001, "DEL") == set()
    _write_all_family_events(tmp_path, {"F1": ["kid", "dad", "mum"], "F2": ["stranger"]},
                             samples, [], nalt, nalt.astype(np.int16), nalt,
                             lambda *a: "SHANK3", impact, _FamilySpill(tmp_path / "inh"),
                             subthreshold=cands)
    page = (tmp_path / "families" / "F1" / "events.html").read_text()
    rows = [r for r in page.split("<tr ")[1:] if 'data-site="sub_' in r]
    assert len(rows) == 1
    row = rows[0][:rows[0].index("</tr>")]
    assert 'data-site="sub_kid_chr22_50731000"' in row and 'data-sub="1"' in row
    assert 'data-pass="0"' in row and "SOUS-SEUIL NDD" in row and "SQ 0.38" in row
    assert "<td>50,731,001</td><td>50,740,000</td>" in row
    # visible by default: the sub-threshold box is ticked next to "PASS uniquement"
    assert 'id="sub" checked' in page


def test_subthreshold_rows_of_a_noisy_genome_are_masked():
    """dataset_paper: 71 sub-threshold NDD rows in one genome, 53 in another —
    noise, not events. Past SUB_MAX_PER_SAMPLE a genome's rows are masked and
    counted for its QC alert; a genome with one row keeps it."""
    from vedette.report import SUB_MAX_PER_SAMPLE, SubCalls, subthreshold_events

    sc = SubCalls()
    n = SUB_MAX_PER_SAMPLE + 1
    sc.add("noisy", [("chr1", "DEL", 1_000_000 * (k + 1), 1_000_000 * (k + 1) + 9_000, 1, 99, 0.1)
                     for k in range(n)])
    sc.add("clean", [("chr2", "DEL", 5_000_000, 5_009_000, 1, 99, 0.4)])
    sc.freeze()

    def impact(*a):
        return 0.9, "+0.90 — élevé · DEL exonique NDD1"

    rows, masked = subthreshold_events(["noisy", "clean"], sc, [], None, impact)
    assert masked == {"noisy": n}
    assert [r["carriers"] for r in rows] == [["clean"]]


def test_large_subthreshold_dup_pieces_are_joined_into_a_row():
    """The v5 IGV review: a de novo 15q11.2q13 duplication at CN 4 was cut
    into LowSQ pieces by the repeats inside the region; joined, it is a row
    whatever its genes. A common CNV it runs through does not make it shared;
    a large DUP found outside the family, a short one with no NDD gene and the
    X of an XXY boy are not rows."""
    from vedette.report import SubCalls, join_dup_pieces, subthreshold_events

    pieces = [("chr15", "DUP", 22_606_000, 23_124_001, 4, 99, 0.02),   # BP1-BP2, common
              ("chr15", "DUP", 23_360_000, 24_427_001, 4, 99, 0.06),
              ("chr15", "DUP", 24_477_000, 30_158_001, 4, 99, 0.03)]
    sc = SubCalls()
    sc.add("kid", pieces + [
        ("chr7", "DUP", 10_000_000, 10_050_001, 3, 99, 0.1),          # 50 kb, no gene
        ("chr8", "DUP", 2_190_000, 2_380_001, 3, 99, 0.15),           # shared outside
        ("chrX", "DUP", 148_173_000, 149_654_001, 3, 99, 0.0),        # XXY: the extra X
    ])
    sc.add("stranger", [("chr15", "DUP", 22_610_000, 23_120_001, 3, 99, 0.1),
                        ("chr8", "DUP", 2_200_000, 2_370_001, 3, 99, 0.2)])
    sc.freeze()
    joined = join_dup_pieces(sc.calls("kid"))
    big = [j for j in joined if j["contig"] == "chr15"]
    assert [(j["start"], j["end"], len(j["pieces"])) for j in big] == [(22_606_000, 30_158_001, 3)]
    assert big[0]["sq"] == pytest.approx(0.06)

    def impact(*a):
        return 0.0, "+0.00 — faible"

    rows, masked = subthreshold_events(["kid", "dad", "mum"], sc, [], None, impact,
                                       sex_aneuploid={"kid": {"X"}})
    assert masked == {}
    assert [(r["contig"], r["pos"], r["end"], r["rescue"]) for r in rows] == [
        ("chr15", 22_606_001, 30_158_001, "SOUS-SEUIL GRANDE DUP")]
    assert "3 appels joints" in rows[0]["inherit"] and rows[0]["gt"][0] == ("0/1", "99", "4")
    # pieces too far apart stay apart: 100 kb gap on a 250 kb span
    two = join_dup_pieces([{"contig": "chr1", "svtype": "DUP", "start": a, "end": b,
                            "cn": 3, "gq": 99, "sq": None}
                           for a, b in ((0, 75_000), (175_000, 250_000))])
    assert len(two) == 2


def test_inherited_subthreshold_dup_is_one_row():
    """PMS v1.6.0: a parent and child carrying the same large sub-threshold
    duplication came out as two rows; it is one row with both carriers."""
    from vedette.report import SubCalls, subthreshold_events

    sc = SubCalls()
    sc.add("kid", [("chr4", "DUP", 77_824_000, 78_061_001, 3, 99, 0.1)])
    sc.add("dad", [("chr4", "DUP", 77_825_000, 78_061_001, 3, 95, 0.3)])
    sc.freeze()
    rows, _ = subthreshold_events(["kid", "dad", "mum"], sc, [], None, lambda *a: (0.0, ""))
    assert len(rows) == 1
    r = rows[0]
    assert r["carriers"] == ["kid", "dad"] and r["sqmax"] == "0.30"
    assert (r["pos"], r["end"]) == (77_824_001, 78_061_001)
    assert [g[0] for g in r["gt"]] == ["0/1", "0/1", "—"]


def test_family_export_holds_rare_non_pass_sites(tmp_path):
    """#10: the page shows PASS sites only (--events-pass-only), the family
    export also the non-PASS sites carried by at most 20 % of the cohort;
    those carried by more (dataset_paper: 167k of 190k non-PASS sites) are
    left out."""
    import csv

    import numpy as np

    from vedette.report import _FamilySpill, _write_all_family_events

    samples = ["kid", "dad", "mum"]

    def site(i, start, filt, cf):
        return {"id": f"s{i}", "contig": "chr1", "start": start, "end": start + 5_000,
                "svtype": "DEL", "size": 5_000, "filter": filt, "rescue": "", "idx": i,
                "sqmax": None, "carrier_frac": cf, "ref_carrier_frac": None}

    sites = [site(0, 10_000, "PASS", 0.05), site(1, 100_000, "CohortArtifact", 0.12),
             site(2, 200_000, "HighCarrierFrac", 0.90), site(3, 300_000, "CommonRef", None)]
    nalt = np.array([[1, 0, 0]] * 4, dtype=np.int8)
    gq = np.full((4, 3), 99, dtype=np.int16)
    cn = np.array([[1, 2, 2]] * 4, dtype=np.int8)

    def run(**kw):
        _write_all_family_events(tmp_path, {"F1": samples}, samples, sites, nalt, gq, cn,
                                 lambda *a: "", lambda *a: (0, ""), _FamilySpill(tmp_path / "inh"),
                                 pass_only=True, **kw)
        with open(tmp_path / "families" / "F1" / "F1.vedette.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        page = (tmp_path / "families" / "F1" / "events.html").read_text()
        return [(r["row_id"], r["filter"], r["kind"]) for r in rows], page

    rows, page = run()
    assert rows == [("s0", "PASS", "site"), ("s1", "CohortArtifact", "site")]
    assert "CohortArtifact" not in page and "100,001" not in page
    rows, _ = run(export_nonpass_max_cf=0)
    assert rows == [("s0", "PASS", "site")]


def test_large_family_gets_no_dead_events_link():
    """A 38-member dataset_paper pedigree had a link to an events page that
    is never written (pages stop at _EVENTS_MAX_MEMBERS members)."""
    from vedette.report import _EVENTS_MAX_MEMBERS, _events_link

    assert 'href="families/F1/events.html"' in _events_link("F1", 3)
    big = _events_link("F2", _EVENTS_MAX_MEMBERS + 18)
    assert "href" not in big and f"{_EVENTS_MAX_MEMBERS + 18} membres" in big


def test_wave_backfill_on_a_store_without_score(tmp_path):
    """Stores called before schema 0.4 carry no wave_score (2,517 of 2,670 on
    dataset_paper): the backfill recomputes it from the stored depth, equal to
    what the call computes, and leaves an existing score alone."""
    from vedette.extract.store import EvidenceStore, EvidenceStoreWriter
    from vedette.qc import backfill_wave, sample_wave_score

    rng = np.random.default_rng(0)
    path = tmp_path / "S.evidence.zarr"
    w = EvidenceStoreWriter(path, sample="S", reference_id="x", bin_size=10_000)
    ratio, mask = {}, {}
    for k in range(1, 13):
        n = 2_000
        x = 1.0 + 0.05 * np.sin(np.arange(n) / 40) + rng.normal(0, 0.02, n)
        ratio[f"chr{k}"], mask[f"chr{k}"] = x, np.zeros(n, dtype=bool)
        w.write_depth(f"chr{k}", (30 * x).astype(np.float32), (30 * x).astype(np.float32), mask[f"chr{k}"])
    expected, _ = sample_wave_score(ratio, mask, 10_000)
    assert "wave_score" not in EvidenceStore(path).meta
    score, written = backfill_wave(path)
    assert written and score == pytest.approx(expected, abs=1e-6)
    meta = EvidenceStore(path).meta
    assert meta["wave_score"] == pytest.approx(expected, abs=1e-6) and len(meta["wave_by_contig"]) == 12
    assert backfill_wave(path) == (pytest.approx(score), False)


def test_denovo_alerts_split_the_reviewed_false_from_the_true():
    """IGV review (2026-10-13): the 7 false high-impact de novo rows were
    177 bp - 3.1 kb, one 3.1 kb site in 3 families; true ones >= 10.5 kb,
    including 22q13 deletions shared by up to 40 PMS families."""
    from vedette.report import denovo_alerts

    for size, n_fam in ((177, 2), (250, 1), (313, 4), (447, 2), (3131, 3)):
        assert denovo_alerts(size, n_fam)
    for size, n_fam in ((10_511, 1), (81_051, 40), (365_000, 30), (623_000, 2), (8_905_000, 2)):
        assert denovo_alerts(size, n_fam) == []
    assert denovo_alerts(3131, 3) == ["de novo récurrent (3 familles)"]


def test_no_subthreshold_rows_for_low_quality_or_wavy_genomes():
    """Sub-threshold rows were true 1/10 on old IP genomes and 1/13 on wavy
    ones (IGV review 2026-10-13): their candidates are skipped."""
    from vedette.report import SubCalls, subthreshold_events

    sc = SubCalls()
    for s in ("old_genome", "wavy", "good"):
        sc.add(s, [("chr2", "DEL", 5_000_000, 5_009_000, 1, 99, 0.4)])
    sc.freeze()

    def impact(*a):
        return 0.9, "+0.90 — élevé · DEL exonique NDD1"

    rows, _ = subthreshold_events(["old_genome", "wavy", "good"], sc, [], None, impact,
                                  skip={"old_genome", "wavy"})
    assert [r["carriers"] for r in rows] == [["good"]]


def test_nested_de_novo_sites_of_one_child_are_one_event():
    """A PMS proband's 8.2 Mb terminal deletion was carried at 61 nested
    cohort sites (one per breakpoint cluster of other probands), each listed
    as a de novo row; the child's own calls chain them into one event."""
    from vedette.report import collapse_denovo

    def row(site, start, end, child="kid"):
        return {"site_id": site, "contig": "chr22", "start": str(start), "end": str(end),
                "svtype": "DEL", "child": child}

    rows = [row("s1", 42_559_000, 46_859_002), row("s2", 48_992_000, 50_800_002),
            row("s3", 50_669_805, 50_750_857), row("s4", 10_000_000, 10_050_000),
            row("s6", 42_761_000, 50_750_857),  # another proband's deletion, spanning the join
            row("s5", 48_992_000, 50_800_002, child="sib")]
    # the child's own PASS calls: five pieces 2 kb apart, 42.56-50.80 Mb
    pieces = [(42_559_001, 46_859_001), (46_861_001, 48_986_001), (48_991_001, 49_384_001),
              (49_386_001, 50_344_001), (50_346_001, 50_799_001)]
    idx = {("chr22", "DEL"): [(a, b, True, "kid", True) for a, b in pieces]}
    evs = collapse_denovo(rows, own_idx_of=lambda c: idx if c == "kid" else {})
    kid = sorted((e for e in evs if e["child"] == "kid"), key=lambda e: e["pos"])
    assert [(e["pos"], e["end"], len(e["sites"])) for e in kid] == [
        (10_000_001, 10_050_000, 1), (42_559_001, 50_799_001, 4)]
    assert [len(e["sites"]) for e in evs if e["child"] == "sib"] == [1]


def _cn_store(run_dir, sample, segments, n_bins=3_000, masked=(), depth=40.0, seed=0):
    """Store with chr22 at CN 2 except `segments` [(start_bin, end_bin, cn)]."""
    from vedette.extract.store import EvidenceStoreWriter

    rng = np.random.default_rng(seed)
    cn = np.full(n_bins, 2.0)
    for a, b, c in segments:
        cn[a:b] = c
    corr = (depth * cn / 2 + rng.normal(0, 1.5, n_bins)).astype(np.float32)
    mask = np.zeros(n_bins, dtype=bool)
    for a, b in masked:
        mask[a:b] = True
    d = run_dir / "samples" / sample
    d.mkdir(parents=True, exist_ok=True)
    w = EvidenceStoreWriter(d / f"{sample}.evidence.zarr", sample=sample, reference_id="x", bin_size=1000)
    w.write_depth("chr22", corr, corr, mask)
    w.set_meta(mean_depth=depth, ploidy={"chr22": 2})


def test_row_bounds_follow_the_carriers_depth(tmp_path):
    """IGV review (2026-10-13): row starts off by 0.3-1.4 Mb — a cohort site
    from another proband's breakpoints, a chain cut at one piece. The child's
    deletion runs 1,000-2,800 kb here, crossing a masked stretch."""
    from vedette.report import _DepthCache, refine_bounds

    _cn_store(tmp_path, "kid", [(1_000, 2_800, 1.0)], masked=[(1_500, 1_560)])
    _cn_store(tmp_path, "dad", [], seed=1)
    _cn_store(tmp_path, "mum", [], seed=2)
    depth = _DepthCache(tmp_path)
    # row starting after the masked stretch, as when the chain was cut there
    pos, end = refine_bounds(depth, "chr22", 1_600_001, 2_800_000, "DEL", ["kid"], ["dad", "mum"])
    assert abs(pos - 1_000_001) <= 3_000 and abs(end - 2_800_000) <= 3_000
    # row starting 300 kb too early (another proband's site)
    pos, _ = refine_bounds(depth, "chr22", 700_001, 2_800_000, "DEL", ["kid"], ["dad", "mum"])
    assert abs(pos - 1_000_001) <= 3_000


def test_row_bounds_do_not_swallow_an_event_the_parent_carries(tmp_path):
    """The child's de novo deletion (1,500-2,000 kb) abuts a deletion it
    inherited with the mother (1,000-1,500 kb): the de novo row stays put."""
    from vedette.report import _DepthCache, refine_bounds

    _cn_store(tmp_path, "kid", [(1_000, 2_000, 1.0)])
    _cn_store(tmp_path, "mum", [(1_000, 1_500, 1.0)], seed=1)
    _cn_store(tmp_path, "dad", [], seed=2)
    depth = _DepthCache(tmp_path)
    assert refine_bounds(depth, "chr22", 1_500_001, 2_000_000, "DEL", ["kid"], ["mum", "dad"]) is None


def test_no_bound_refinement_for_duplications_or_very_wavy_carriers(tmp_path):
    """Second IGV review (2026-10-14): both refined DUPs were wrong (a wave
    read as a gain, 112 kb -> 23 Mb) and refinement was a coin toss on
    carriers with wave >= 0.2."""
    from vedette.report import REFINE_MAX_WAVE, _DepthCache, refine_bounds

    _cn_store(tmp_path, "kid", [(1_000, 2_800, 3.0)])
    _cn_store(tmp_path, "dad", [], seed=1)
    _cn_store(tmp_path, "mum", [], seed=2)
    assert refine_bounds(_DepthCache(tmp_path), "chr22", 1_600_001, 2_800_000, "DUP",
                         ["kid"], ["dad", "mum"]) is None
    import zarr

    _cn_store(tmp_path / "w", "kid", [(1_000, 2_800, 1.0)])
    _cn_store(tmp_path / "w", "dad", [], seed=1)
    zarr.open_group(str(tmp_path / "w" / "samples" / "kid" / "kid.evidence.zarr"), mode="r+").attrs[
        "wave_score"] = REFINE_MAX_WAVE + 0.05
    assert refine_bounds(_DepthCache(tmp_path / "w"), "chr22", 1_600_001, 2_800_000, "DEL",
                         ["kid"], ["dad"]) is None


def test_mosaic_fragments_fold_into_the_mosaic_row(tmp_path):
    """Second IGV review: a clonal genome's arm-level losses were spread over
    CN 1 fragment rows beside the mosaic row. Fragments of the sample alone at
    the mosaic's level fold in; a deeper deletion inside the mosaic, and an
    event the parent also carries, stay rows of their own."""
    from vedette.report import _DepthCache, fold_mosaic_fragments

    _cn_store(tmp_path, "kid", [(1_000, 2_500, 1.33), (1_800, 1_900, 1.0)])
    mosaic = {"own": "mosaïque", "carriers": ["kid"], "contig": "chr22", "svtype": "DEL",
              "pos": 1_000_001, "end": 2_500_000, "mosaic_cn": 1.33, "n_sites": 1,
              "inherit": "mosaïque 67 % (CN 1.33)", "id": "mosaic_kid"}

    def row(i, a, b, carriers=("kid",)):
        return {"own": "PASS", "carriers": list(carriers), "contig": "chr22", "svtype": "DEL",
                "pos": a, "end": b, "id": i}

    rows = [mosaic, row("frag1", 1_100_001, 1_300_000), row("frag2", 1_400_001, 1_700_000),
            row("deep", 1_800_001, 1_900_000), row("shared", 2_000_001, 2_200_000, ("kid", "mum")),
            row("outside", 2_600_001, 2_700_000)]
    out = fold_mosaic_fragments(rows, _DepthCache(tmp_path))
    assert [e["id"] for e in out] == ["mosaic_kid", "deep", "shared", "outside"]
    assert out[0]["n_sites"] == 3 and "2 fragments regroupés" in out[0]["inherit"]
    assert {"frag1", "frag2"} <= set(out[0]["_ids"])


def test_large_duplication_keeps_pass_when_the_allelic_balance_shows_a_gain():
    """A trisomy 21, a 15q11-q13 duplication, a 7q11.23 duplication and a
    trisomy 13 were all LowSQ (SQ 0.01-0.09): the model never saw a real large
    depth-only gain. Their heterozygous sites sit at 1/3-2/3 (mean |VAF - 0.5|
    ~0.14); a coverage wave leaves them at 0.5, and segmental duplications
    fake the imbalance, so those stay rejected."""
    from vedette.config import load_config
    from vedette.filter.pipeline import sq_override

    f = load_config(None).filter
    gain = {"cn": 3, "baf_n": 260, "baf_dev": 0.142, "segdup_ovl": 0.0, "maskfrac": 0.0,
            "fc_chrom": 1.4}
    assert sq_override("DUP", 994_000, gain, f) == "LARGE_DUP_BAF"
    assert sq_override("DUP", 994_000, dict(gain, baf_dev=0.06), f) is None   # a wave
    assert sq_override("DUP", 994_000, dict(gain, segdup_ovl=0.8), f) is None  # paralogs
    assert sq_override("DUP", 994_000, dict(gain, baf_n=12), f) is None        # too few hets
    assert sq_override("DUP", 60_000, gain, f) is None                         # model's range
    assert sq_override("DUP", 994_000, dict(gain, cn=2), f) is None
    assert sq_override("DEL", 6_427_000, dict(gain, fc_chrom=0.528, cn=1), f) == "LARGE_DEL_DEPTH"


def test_constitutional_segments_are_rows_unless_calls_cover_them():
    """A trisomy 21 came back as LowSQ pieces and was on no family page: the
    scan's whole-copy segment becomes a row. A PMS terminal deletion already
    shown by PASS rows keeps those rows (finer bounds than 1 Mb windows)."""
    from vedette.report import constitutional_events

    seg = {"contig": "chr21", "start": 15_000_000, "end": 46_000_000, "svtype": "DUP",
           "cn": 3.07, "fraction": 1.0, "mosaic": False}
    deletion = {"contig": "chr22", "start": 46_000_000, "end": 50_818_468, "svtype": "DEL",
                "cn": 1.06, "fraction": 0.94, "mosaic": False}
    pass_row = {"contig": "chr22", "svtype": "DEL", "pos": 46_100_001, "end": 50_800_000,
                "filter": "PASS", "carriers": ["kid"], "own": "PASS"}
    rows = constitutional_events(["kid", "mum"], {"kid": [seg, deletion]}, [pass_row])
    assert [r["contig"] for r in rows] == ["chr21"]
    r = rows[0]
    assert r["rescue"] == "CHROMOSOMIQUE" and r["own"] == "chromosomique"
    assert r["filter"] == "PASS" and "gain constitutionnel (CN 3.07" in r["inherit"]
    assert (r["pos"], r["end"]) == (15_000_001, 46_000_000)


def test_missing_evidence_stores_warn_then_fail(tmp_path, caplog):
    """A re-filter run (--calls_dir) published no store under samples/<S>/:
    the report dropped its karyotype section and wave alerts silently
    (2670 samples, 2026-11-02). Some missing: a warning naming the count and
    the lost sections; none at all: an error before the VCF is decoded."""
    import logging

    import pysam

    from vedette.report import check_evidence_stores, generate_report

    samples = ["a", "b", "c"]
    (tmp_path / "samples" / "a" / "a.evidence.zarr").mkdir(parents=True)
    with caplog.at_level(logging.WARNING, logger="vedette.report"):
        assert check_evidence_stores(tmp_path, samples) == ["b", "c"]
    msg = caplog.text
    assert "2/3" in msg and "Chromosomes sexuels" in msg and "vague" in msg and "b, c" in msg

    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="vedette.report"):
        assert check_evidence_stores(tmp_path, ["a"]) == []
    assert not caplog.records

    with pytest.raises(FileNotFoundError, match="aucun store"):
        check_evidence_stores(tmp_path, ["b", "c"])
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="vedette.report"):
        assert check_evidence_stores(tmp_path, ["b", "c"], require=False) == ["b", "c"]
    assert "2/2" in caplog.text

    # generate_report checks the header's samples before anything else
    run = tmp_path / "run"
    (run / "cohort").mkdir(parents=True)
    hdr = pysam.VariantHeader()
    hdr.add_line("##source=vedette-1.3.3")  # legacy: POS = start + 1, END = end + 1
    hdr.contigs.add("chr1", length=1_000_000)
    hdr.add_sample("kid")
    with pysam.VariantFile(str(run / "cohort" / "cohort.vcf.gz"), "wz", header=hdr):
        pass
    with pytest.raises(FileNotFoundError, match="--allow-missing-stores"):
        generate_report(run, None)


def test_inherit_label_karyotype():
    from vedette.report import _inherit_fr

    assert _inherit_fr("denovo") == "de novo"
    assert _inherit_fr("karyotype_XXY") == "caryotype XXY"
    assert "de novo" not in _inherit_fr("karyotype_XXY")   # no de novo alert


def test_child_quality_names_what_weakens_a_genome():
    from vedette.report import WAVE_WARN, child_quality

    waves = {"A": WAVE_WARN + 0.01, "B": 0.001, "IP00X": None}
    assert child_quality("A", waves) == f"wavy {WAVE_WARN + 0.01:.3f}"
    assert child_quality("B", waves) == ""
    assert child_quality("IP00X", waves, ("IP",), {"IP00X": 6.2}) == "génome ancien, charge 6.2×"
