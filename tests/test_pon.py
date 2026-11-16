import numpy as np
import pysam
import pytest

from svengers.cohort.pon import PanelOfNormals, artifacts_bed, build_pon
from svengers.extract.store import EvidenceStoreWriter

N_BINS = 1000
MEAN_DEPTH = 30.0


def _make_store(tmp_path, sample, depth_ratio, sex="female", contig="chr1"):
    """depth_ratio: per-bin copy ratio; store holds corr = ratio * mean_depth."""
    path = tmp_path / f"{sample}.zarr"
    w = EvidenceStoreWriter(path, sample=sample, reference_id="testref", bin_size=1000)
    ploidy = {"chr1": 2, "chrX": 1 if sex == "male" else 2}
    w.set_meta(mean_depth=MEAN_DEPTH, coverage=MEAN_DEPTH, ploidy=ploidy)
    corr = (np.asarray(depth_ratio, dtype=np.float32) * MEAN_DEPTH)
    w.write_depth(contig, corr, corr, np.zeros(len(corr), dtype=bool))
    return path


@pytest.fixture()
def cohort_stores(tmp_path):
    rng = np.random.default_rng(0)
    paths = []
    for i in range(6):
        ratio = rng.normal(1.0, 0.05, N_BINS)
        ratio[100] = rng.normal(0.5, 0.02)      # systematic reference bias
        ratio[200] = 0.2 if i % 2 else 2.5      # population-bimodal bin
        if i == 0:  # sample 0 carries a real het DEL at bins 300-320
            ratio[300:320] = rng.normal(0.5, 0.05, 20)
        paths.append(_make_store(tmp_path, f"S{i}", ratio))
    return paths


def test_build_and_normalize(tmp_path, cohort_stores):
    pon_path = build_pon(cohort_stores, tmp_path / "pon.zarr", min_samples=5)
    pon = PanelOfNormals(pon_path)
    median, mad = pon.profile("chr1")
    assert median[100] == pytest.approx(0.5, abs=0.1)     # bias captured
    assert median[500] == pytest.approx(1.0, abs=0.05)
    assert mad[200] > 0.3                                  # variability captured
    # the real DEL in 1/6 samples must NOT shift the median
    assert np.median(median[300:320]) == pytest.approx(1.0, abs=0.08)

    # normalize the carrier sample
    from svengers.extract.store import EvidenceStore

    s0 = EvidenceStore(cohort_stores[0])
    _raw, corr, mask = s0.depth("chr1")
    ratio = corr / MEAN_DEPTH
    norm, new_mask = pon.normalize("chr1", ratio, mask)
    assert norm[100] == pytest.approx(1.0, abs=0.2)        # bias divided out
    assert new_mask[200] and not new_mask[500]             # variable bin masked
    assert np.median(norm[300:320]) == pytest.approx(0.5, abs=0.1)  # DEL preserved


def test_sex_stratified_chrx(tmp_path):
    rng = np.random.default_rng(1)
    paths = []
    for i in range(4):
        r = rng.normal(0.5, 0.03, N_BINS)  # male X ratio
        paths.append(_make_store(tmp_path, f"M{i}", r, sex="male", contig="chrX"))
    for i in range(4):
        r = rng.normal(1.0, 0.03, N_BINS)
        paths.append(_make_store(tmp_path, f"F{i}", r, sex="female", contig="chrX"))
    pon = PanelOfNormals(build_pon(paths, tmp_path / "pon.zarr", min_samples=8))
    med_m, _ = pon.profile("chrX", "male")
    med_f, _ = pon.profile("chrX", "female")
    assert np.median(med_m) == pytest.approx(0.5, abs=0.05)
    assert np.median(med_f) == pytest.approx(1.0, abs=0.05)


def test_male_chrx_stays_at_its_baseline(tmp_path):
    """A male chrX divided by the male panel must come back at ratio 0.5 (CN1),
    not 1.0: every dataset_paper male had his whole chrX called a DUP."""
    rng = np.random.default_rng(2)
    paths = []
    for sex, base in (("male", 0.5), ("female", 1.0)):
        for i in range(4):
            r = rng.normal(base, 0.03, N_BINS)
            r[:50] = rng.normal(1.0, 0.03, 50)          # pseudo-autosomal: 2 copies in both
            paths.append(_make_store(tmp_path, f"{sex}{i}", r, sex=sex, contig="chrX"))
    pon = PanelOfNormals(build_pon(paths, tmp_path / "pon.zarr", min_samples=8))
    mask = np.zeros(N_BINS, dtype=bool)
    male = rng.normal(0.5, 0.03, N_BINS)
    male[:50] = 1.0
    male[600:650] = 0.0                                  # hemizygous deletion
    norm, _ = pon.normalize("chrX", male, mask, sex="male", level=0.5)
    assert np.median(norm[100:600]) == pytest.approx(0.5, abs=0.05)
    assert np.median(norm[:50]) == pytest.approx(0.5, abs=0.05)
    assert np.max(norm[600:650]) == pytest.approx(0.0, abs=0.01)
    female = rng.normal(1.0, 0.03, N_BINS)
    norm, _ = pon.normalize("chrX", female, mask, sex="female", level=1.0)
    assert np.median(norm) == pytest.approx(1.0, abs=0.05)
    # chrY of a female (ploidy 0): left as it is
    norm, new_mask = pon.normalize("chrX", female, mask, sex="female", level=0.0)
    assert np.array_equal(norm, female) and not new_mask.any()


def test_bin_size_guard(tmp_path, cohort_stores):
    pon = PanelOfNormals(build_pon(cohort_stores, tmp_path / "pon.zarr", min_samples=5))
    assert pon.bin_size == 1000


def test_svd_denoising_removes_batch_wave_keeps_cnv(tmp_path):
    """A shared low-rank batch effect must be removed; a sample-private DEL kept."""
    rng = np.random.default_rng(3)
    bins = np.arange(N_BINS)
    wave = 0.25 * np.sin(2 * np.pi * bins / 200)  # systematic pattern (log2 units)
    paths = []
    amplitudes = rng.normal(0.0, 1.0, 30)
    for i in range(30):
        ratio = rng.normal(1.0, 0.03, N_BINS) * np.exp2(amplitudes[i] * wave)
        paths.append(_make_store(tmp_path, f"S{i}", ratio))

    # target: extreme batch amplitude + a real het DEL at bins 400-450
    target_ratio = rng.normal(1.0, 0.03, N_BINS) * np.exp2(2.5 * wave)
    target_ratio[400:450] *= 0.5
    target_path = _make_store(tmp_path, "TARGET", target_ratio)

    pon = PanelOfNormals(
        build_pon(paths, tmp_path / "pon.zarr", min_samples=8, denoise_k=3)
    )
    from svengers.extract.store import EvidenceStore

    s = EvidenceStore(target_path)
    _raw, corr, mask = s.depth("chr1")
    ratio = corr / MEAN_DEPTH

    raw_out, _ = pon.normalize("chr1", ratio, mask, denoise=False)
    svd_out, _ = pon.normalize("chr1", ratio, mask, denoise=True)

    clean = np.ones(N_BINS, dtype=bool)
    clean[400:450] = False
    # residual wave clearly visible without SVD, mostly gone with it
    resid_raw = np.std(np.log2(raw_out[clean]))
    resid_svd = np.std(np.log2(svd_out[clean]))
    assert resid_svd < 0.5 * resid_raw
    # the DEL survives denoising
    assert np.median(svd_out[400:450]) == pytest.approx(0.5, abs=0.08)


def test_artifacts_bed(tmp_path):
    h = pysam.VariantHeader()
    h.add_line("##source=svengers-1.3.3")  # legacy: POS = start + 1, END = end + 1
    h.contigs.add("chr1", length=1_000_000)
    for line in (
        ("SVTYPE", "1", "String", "t"), ("END", "1", "Integer", "e"),
        ("AF", "A", "Float", "af"), ("HWE_P", "1", "Float", "p"),
        ("GTR_DISC", "1", "Float", "d"), ("NSAMPLES", "1", "Integer", "n"),
    ):
        h.info.add(*line)
    vcf_path = tmp_path / "cohort.vcf.gz"
    with pysam.VariantFile(str(vcf_path), "wz", header=h) as out:
        specs = [
            (1000, 2000, 0.10, 1e-9, 0.0, 10),   # artifact: common + HWE failure
            (5000, 6000, 0.10, 0.5, 0.0, 10),    # common but HWE fine -> keep
            (8000, 9000, 0.01, 0.5, 0.8, 5),     # discordant nominators -> artifact
            (12000, 13000, 0.01, 0.5, 0.0, 2),   # rare and clean -> keep
        ]
        for start, end, af, hwe, disc, ns in specs:
            rec = out.new_record(
                contig="chr1", start=start, stop=end, alleles=("N", "<DEL>")
            )
            rec.info.update(
                SVTYPE="DEL", AF=(af,), HWE_P=hwe, GTR_DISC=disc, NSAMPLES=ns
            )
            out.write(rec)
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)

    out_bed = tmp_path / "artifacts.bed"
    n = artifacts_bed(vcf_path, out_bed)
    lines = out_bed.read_text().strip().splitlines()
    assert n == 2 and len(lines) == 2
    starts = {int(line.split("\t")[1]) for line in lines}
    assert starts == {1000, 8000}
