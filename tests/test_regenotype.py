"""Re-genotyping a sample's calls from its store after the BAF is rewritten."""

import numpy as np
import pysam

from vedette.extract.store import EvidenceStoreWriter
from vedette.genotype.regenotype import regenotype_calls


def _store(path, vaf):
    w = EvidenceStoreWriter(path, sample="S1", reference_id="x", bin_size=1000)
    w.set_meta(mean_depth=30.0, coverage=30.0, ploidy={"chr1": 2})
    depth = np.full(2000, 30.0, np.float32)
    depth[500:600] = 15.0                                   # a heterozygous deletion
    w.write_depth("chr1", depth, depth, np.zeros(2000, bool))
    pos = np.arange(100_000, 1_000_000, 5_000, dtype=np.uint32)
    w.write_baf("chr1", pos, vaf(pos).astype(np.float32), np.full(len(pos), 30, np.uint16))


def _calls(path):
    hdr = pysam.VariantHeader()
    hdr.add_line("##vedette_coords=padded")
    hdr.contigs.add("chr1", length=2_000_000)
    hdr.add_line('##ALT=<ID=DEL,Description="Deletion">')
    for k, n, t in (("SVTYPE", "1", "String"), ("END", "1", "Integer"), ("SOURCE", "1", "String")):
        hdr.info.add(k, n, t, k)
    for k, t in (("GT", "String"), ("GQ", "Integer"), ("CN", "Integer"), ("CNQ", "Float"),
                 ("RD", "Float"), ("PE", "Integer"), ("SR", "Integer"), ("CG", "Integer"),
                 ("BAFN", "Integer"), ("BAFDEV", "Float")):
        hdr.formats.add(k, "1", t, k)
    hdr.add_sample("S1")
    with pysam.VariantFile(str(path), "wz", header=hdr) as vf:
        rec = vf.new_record(contig="chr1", start=499_999, stop=600_000, alleles=("N", "<DEL>"),
                            id="vedette_DEL_1")
        rec.info["SVTYPE"], rec.info["SOURCE"] = "DEL", "DEPTH"
        s = rec.samples["S1"]
        s["GT"], s["GQ"], s["PE"], s["SR"], s["CG"] = (0, 1), 7, 0, 0, 0
        s["BAFN"], s["BAFDEV"] = 2, 0.05                    # stale: an old, sparse BAF
        vf.write(rec)


def test_regenotype_takes_the_rewritten_baf(tmp_path):
    """A heterozygous deletion has no heterozygous site: with a BAF that has
    them elsewhere and none inside, BAFN is 0 there; GT/GQ/CN come from depth
    and BAF as at calling time, the junction support is kept."""
    store = tmp_path / "S1.evidence.zarr"
    # hets outside the deletion, hom-alt inside (loss of heterozygosity)
    _store(store, lambda pos: np.where((pos >= 500_000) & (pos < 600_000), 1.0, 0.5))
    vcf, out = tmp_path / "S1.vedette.vcf.gz", tmp_path / "S1.regeno.vcf.gz"
    _calls(vcf)
    c = regenotype_calls(vcf, store, out)
    assert (c["calls"], c["regenotyped"], c["unchanged"]) == (1, 1, 0)
    rec = next(iter(pysam.VariantFile(str(out))))
    s = rec.samples[0]
    assert s["BAFN"] == 0 and s.get("BAFDEV") is None
    assert s["CN"] == 1 and tuple(s["GT"]) == (0, 1) and s["GQ"] > 7
    assert abs(s["RD"] - 0.5) < 0.01 and s["PE"] == 0
    assert "##vedette_regenotyped=" in str(pysam.VariantFile(str(out)).header)


def test_cohort_genotype_reads_the_original_baf(tmp_path):
    """`vedette baf` keeps the store's original BAF as baf_track. The cohort
    genotyper reads it unless genotype.rewritten_baf is set (the rewritten
    BAF has spurious hets in paralogous loci, each one a mismap penalty under
    a deletion); regenotype reads the rewritten one, its purpose."""
    import dataclasses

    import zarr

    from vedette.cohort.genotyper import StoreGenotyper
    from vedette.config import Config
    from vedette.extract.store import EvidenceStore, write_baf_group

    path = tmp_path / "S1.evidence.zarr"
    _store(path, lambda pos: np.where((pos >= 500_000) & (pos < 600_000), 1.0, 0.5))
    (path / "chr1" / "baf").rename(path / "chr1" / "baf_track")
    pos = np.arange(100_000, 1_000_000, 5_000, dtype=np.uint32)          # hets everywhere
    write_baf_group(zarr.open_group(str(path), mode="r+")["chr1"], pos,
                    np.full(len(pos), 0.5, np.float32), np.full(len(pos), 30, np.uint16))
    store = EvidenceStore(path)

    def baf_n(config, **kw):
        return StoreGenotyper(store, config, **kw).genotype_site(
            "s1", "chr1", 500_000, 600_000, "DEL")["baf_n"]

    cfg = Config()
    assert baf_n(cfg) == 0
    rewritten = dataclasses.replace(cfg, genotype=dataclasses.replace(cfg.genotype, rewritten_baf=True))
    assert baf_n(rewritten) == 20
    assert baf_n(cfg, rewritten_baf=True) == 20
    # a store never rewritten has no baf_track: the same BAF either way
    plain = tmp_path / "S2.evidence.zarr"
    _store(plain, lambda pos: np.full(len(pos), 0.5))
    assert StoreGenotyper(EvidenceStore(plain), cfg).genotype_site(
        "s1", "chr1", 500_000, 600_000, "DEL")["baf_n"] == 20
