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
    pts = lambda a, b, t="DEL": idx.impact("chr22", a, b, t)[0]
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


def test_denovo_outlier_guard():
    from vedette.report import denovo_outliers

    dn = {
        "F1": [{"child": "kidA"}] * 3,
        "F2": [{"child": "kidB"}] * 25 + [{"child": "kidC"}],
        "F3": [{"child": "kidB"}] * 5,   # counted across families
    }
    out = denovo_outliers(dn, threshold=20)
    assert out == {"kidB": 30}
    assert denovo_outliers({}, 20) == {}


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
    assert by_start[1_000_000]["end"] == 1_020_001  # inside: untouched


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
    assert "<td>45,146,001</td><td>48,986,001</td>" in page
    assert "44,116,001" not in page
    assert (45_146_000, 48_986_001) in seen      # genes of the own span


def test_carrier_requires_sufficient_gq():
    """Re-genotyping emits GT>0 with a low GQ far too readily: on PMS the
    apparent de novo rate — a genotype-error proxy — is 16.7% for DEL
    carriers below GQ 20 against 0.8% above (13.9% vs 1.9% for DUP), so a
    low-GQ alternative genotype must not make the sample a carrier."""
    import numpy as np

    from vedette.report import _write_all_family_events, _FamilySpill

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
    import pysam

    from vedette.config import load_config
    from vedette.filter.pipeline import run_filter

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
