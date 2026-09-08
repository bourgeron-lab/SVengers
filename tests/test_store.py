import numpy as np
import pytest

from vedette.extract.store import EvidenceStore, EvidenceStoreWriter, PE_DEL, PE_DUP


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
