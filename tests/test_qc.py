import numpy as np

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


def test_event_pages_by_chr(tmp_path):
    from vedette.report import _write_event_pages

    members = ["C000EZH", "C000EZI", "C000EZJ"]
    ev = {
        "contig": "chr1", "pos": 1000, "end": 2000, "svtype": "DEL", "size": 1000,
        "filter": "PASS", "rescue": "", "genes": "GENE1", "impact": "3 — exonique, NDD-candidat",
        "carriers": ["C000EZJ"], "gt": [("0/0", 99, "9.0"), ("0/0", 99, ""), ("0/1", 45, "7.2")],
    }
    ev2 = dict(ev, contig="chr2", svtype="DUP", filter="LowSQ", carriers=[], rescue="RECURRENT")
    _write_event_pages(tmp_path / "events", "Famille test", members, [ev, ev2])

    idx = (tmp_path / "events" / "index.html").read_text()
    assert "1 DEL + 1 DUP" in idx and 'href="chr1.html"' in idx and 'href="chr2.html"' in idx
    chr1 = (tmp_path / "events" / "chr1.html").read_text()
    assert "GENE1" in chr1 and "0/1" in chr1 and "NDD-candidat" in chr1
    chr2 = (tmp_path / "events" / "chr2.html").read_text()
    assert "nopass" in chr2 and "RECURRENT" in chr2
