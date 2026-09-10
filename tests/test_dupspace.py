import numpy as np
import pysam
import pytest

from vedette.bench.dupspace import is_tandem_ins


@pytest.fixture(scope="module")
def fasta(tmp_path_factory):
    rng = np.random.default_rng(42)
    seq = "".join(rng.choice(list("ACGT"), size=20_000))
    path = tmp_path_factory.mktemp("ref") / "ref.fa"
    path.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def test_perfect_tandem_dup_detected(fasta):
    pos = 5_000
    dup = fasta.fetch("chr1", pos, pos + 200)  # S duplicates the downstream flank
    assert is_tandem_ins(dup, "chr1", pos, fasta)


def test_upstream_tandem_detected(fasta):
    pos = 8_000
    dup = fasta.fetch("chr1", pos - 150, pos)  # S duplicates the upstream flank
    assert is_tandem_ins(dup, "chr1", pos, fasta)


def test_jittered_tandem_detected(fasta):
    pos = 11_000
    s = list(fasta.fetch("chr1", pos, pos + 120))
    for i in (10, 60, 100):  # a few mismatches
        s[i] = "A" if s[i] != "A" else "C"
    assert is_tandem_ins("".join(s), "chr1", pos, fasta)


def test_novel_sequence_rejected(fasta):
    rng = np.random.default_rng(7)
    novel = "".join(rng.choice(list("ACGT"), size=300))  # ALU-like random insert
    assert not is_tandem_ins(novel, "chr1", 14_000, fasta)


def test_two_copy_tandem_detected(fasta):
    pos = 17_000
    unit = fasta.fetch("chr1", pos, pos + 80)
    assert is_tandem_ins(unit + unit, "chr1", pos, fasta)  # 2x unit insertion
