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
    # every impact label must start with its integer score (the de novo
    # tables sort by that prefix; a bare label crashed the cohort report)
    import pysam

    from vedette.report import _GeneIndex

    gff = tmp_path / "genes.gff"
    lines = [
        'chr1\thavana\tgene\t10000\t20000\t.\t+\t.\tgene_id="ENSG1";gene_name="GENE1";gene_biotype="protein_coding";',
        'chr1\thavana\texon\t10000\t11000\t.\t+\t.\tgene_id="ENSG1";gene_name="GENE1";gene_biotype="protein_coding";',
    ]
    gff.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(gff), preset="gff", force=True)  # bgzips to .gff.gz
    idx = _GeneIndex(str(gff) + ".gz")
    for span in ((10_500, 10_600), (15_000, 16_000), (50_000, 60_000)):
        score, label = idx.impact("chr1", *span)
        assert int(label.split(" ")[0]) == score


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
        "refcf": "0.004", "lcr": "0.00", "segdup": "1.00",
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
