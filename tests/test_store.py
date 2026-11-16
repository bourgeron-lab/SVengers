import numpy as np
import pytest

from svengers.extract.store import PE_DEL, PE_DUP, EvidenceStore, EvidenceStoreWriter


@pytest.fixture()
def store(tmp_path):
    path = tmp_path / "S1.evidence.zarr"
    w = EvidenceStoreWriter(path, sample="S1", reference_id="deadbeef", bin_size=100)
    w.set_meta(mean_depth=30.0, insert_mu=350.0, insert_sigma=70.0, coverage=30.0)
    raw = np.arange(100, dtype=np.float32)
    corr = raw * 1.1
    mask = np.zeros(100, dtype=bool)
    mask[10:20] = True
    w.write_depth("chr1", raw, corr, mask)
    w.write_baf(
        "chr1",
        np.array([500, 1500, 2500], np.uint32),
        np.array([0.5, 0.33, 1.0], np.float32),
        np.array([30, 25, 40], np.uint16),
    )
    w.write_pe(
        "chr1",
        pos1=np.array([1000, 1050, 5000], np.uint32),
        pos2=np.array([4000, 4020, 9000], np.uint32),
        mate_tid=np.array([0, 0, 0], np.int32),
        orient=np.array([PE_DEL, PE_DEL, PE_DUP], np.uint8),
        mapq=np.array([60, 60, 40], np.uint8),
    )
    w.write_sr(
        "chr1",
        pos=np.array([1200, 1201, 4000], np.uint32),
        side=np.array([1, 1, 0], np.uint8),
        mapq=np.array([60, 55, 60], np.uint8),
        cliplen=np.array([30, 25, 40], np.uint16),
        seqs=[b"ACGT" * 8, b"TTTTT" * 5, b"GGGGCCCC" * 5],
    )
    return EvidenceStore(path)


def test_meta(store):
    assert store.sample == "S1"
    assert store.bin_size == 100
    assert store.insert_model().mu == 350.0


def test_depth_slicing(store):
    raw, corr, mask = store.depth("chr1", 1000, 2100)
    assert len(raw) == 11
    assert raw[0] == 10
    assert mask[0] and not mask[-1]  # bins 10-19 masked, bin 20 not
    full_raw, _, _ = store.depth("chr1")
    assert len(full_raw) == 100


def test_pe_range_query(store):
    pe = store.pe("chr1", 900, 1100)
    assert len(pe["pos1"]) == 2
    assert (pe["orient"] == PE_DEL).all()
    assert len(store.pe("chr1", 6000, 7000)["pos1"]) == 0


def test_sr_and_baf_queries(store):
    sr = store.sr("chr1", 1190, 1210)
    assert len(sr["pos"]) == 2
    baf = store.baf("chr1", 0, 2000)
    assert len(baf["pos"]) == 2
    assert baf["vaf"][0] == pytest.approx(0.5)


def test_sr_seq_roundtrip(store):
    sr = store.sr("chr1", 1190, 1210, with_seq=True)
    assert sr["seq"] == [b"ACGT" * 8, b"TTTTT" * 5]  # sorted with their positions
    sr2 = store.sr("chr1", 3900, 4100, with_seq=True)
    assert sr2["seq"] == [b"GGGGCCCC" * 5]


def test_missing_groups_return_empty(store):
    assert store.has_junction_evidence("chr1")
    assert len(store.baf("chr1", 90_000, 99_000)["pos"]) == 0


def test_junction_flag_false_without_pe(tmp_path):
    w = EvidenceStoreWriter(tmp_path / "S2.zarr", "S2", "ref", 1000)
    w.set_meta(mean_depth=30.0)
    w.write_depth("chr1", np.ones(10, np.float32), np.ones(10, np.float32), np.zeros(10, bool))
    s = EvidenceStore(tmp_path / "S2.zarr")
    assert not s.has_junction_evidence("chr1")


def _snv_vcf(path, records):
    """Single-sample bgzipped, indexed VCF: (contig, pos1, ref, alt, qual, filter, ad[, dp])."""
    import pysam

    hdr = pysam.VariantHeader()
    for c in ("chr1", "chr2"):
        hdr.contigs.add(c, length=1_000_000)
    hdr.add_line('##FILTER=<ID=RefCall,Description="reference call">')
    hdr.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    hdr.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">')
    hdr.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    hdr.add_sample("S1")
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        for contig, pos, ref, alt, qual, filt, ad, *dp in records:
            rec = vf.new_record(contig=contig, start=pos - 1, alleles=(ref, alt), qual=qual)
            rec.filter.add(filt)
            rec.samples["S1"]["AD"] = ad
            if dp:
                rec.samples["S1"]["DP"] = dp[0]
            vf.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_baf_backfill_from_a_deepvariant_vcf(store, tmp_path):
    """DeepVariant 1.9 gives many true heterozygotes a QUAL below 30, which the
    GHFC track drops (an old 75x genome kept 1.6k of its 32k chr20 hets): the
    backfill keeps every PASS SNV with DP >= 10, whatever its QUAL, and drops
    RefCall records, indels and low-depth sites. The previous arrays stay as
    baf_track/; a VCF matching none of the store's contigs changes nothing."""
    from svengers.extract.baf import backfill_baf
    from svengers.extract.store import SCHEMA_VERSION

    vcf = tmp_path / "S1.vcf.gz"
    _snv_vcf(vcf, [
        ("chr1", 101, "A", "G", 45.0, "PASS", (14, 16)),     # het, QUAL >= 30
        ("chr1", 201, "C", "T", 12.0, "PASS", (11, 9)),      # het, QUAL < 30: kept
        ("chr1", 301, "G", "A", 0.0, "RefCall", (30, 6)),    # reference call
        ("chr1", 401, "T", "C", 35.0, "PASS", (3, 4)),       # depth 7 < 10
        ("chr1", 451, "T", "C", 35.0, "PASS", (5, 4), 12),   # AD 9, DP 12: kept
        ("chr1", 501, "A", "AT", 40.0, "PASS", (10, 10)),    # indel
        ("chr1", 601, "G", "C", 50.0, "PASS", (0, 25)),      # hom-alt
    ])
    n_sites, n_het, written = backfill_baf(store.path, str(vcf))
    assert (n_sites, n_het, written) == (4, 3, True)
    s = type(store)(store.path)
    b = s.baf("chr1", 0, 10_000)
    assert b["pos"].tolist() == [100, 200, 450, 600]
    assert b["vaf"].tolist() == pytest.approx([16 / 30, 9 / 20, 4 / 9, 1.0])
    assert b["dp"].tolist() == [30, 20, 12, 25]
    assert s.meta["schema_version"] == SCHEMA_VERSION and s.meta["baf_source"].startswith("vcf:")
    import zarr
    old = zarr.open_group(str(store.path / "chr1" / "baf_track"), mode="r")
    assert old["pos"][:].tolist() == [500, 1500, 2500]
    # done once; --force rewrites without touching baf_track
    assert backfill_baf(store.path, str(vcf)) == (0, 0, False)
    assert backfill_baf(store.path, str(vcf), force=True)[2]
    assert zarr.open_group(str(store.path / "chr1" / "baf_track"), mode="r")["pos"][:].tolist() == [500, 1500, 2500]

    other = tmp_path / "other.vcf.gz"
    _snv_vcf(other, [("chr2", 101, "A", "G", 45.0, "PASS", (14, 16))])
    with pytest.raises(ValueError):
        backfill_baf(store.path, str(other), force=True)
    assert type(store)(store.path).baf("chr1", 0, 10_000)["pos"].tolist() == [100, 200, 450, 600]
