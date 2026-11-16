import pytest

from svengers.ped import parse_ped

HEADERLESS = """\
FAM1\tCHILD1\tDAD1\tMOM1\t1\t2
FAM1\tDAD1\t0\tMOM0\t1\t1
FAM1\tMOM1\t0\t0\t2\t1
"""

HEADERED = """\
FID\tsample_id\tFatherBarcode\tMotherBarcode\tSex\tPheno_ASD_curated
220151\tB00ID5W\tC0005J0\tC0005IO\t2\t1
220151\tC0005IO\t0\t0\t2\t-9
220151\tC0005J0\t0\t0\t1\t-9
"""


def _write(tmp_path, content, name="ped.tsv"):
    p = tmp_path / name
    p.write_text(content)
    return p


def test_headerless(tmp_path):
    ped = parse_ped(_write(tmp_path, HEADERLESS))
    assert len(ped.individuals) == 3
    assert ped.individuals["CHILD1"].father_id == "DAD1"
    assert ped.individuals["CHILD1"].sex == 1
    # MOM0 referenced but absent -> DAD1 is still a founder
    founders = {i.sample_id for i in ped.founders()}
    assert founders == {"DAD1", "MOM1"}
    trios = ped.trios()
    assert len(trios) == 1
    assert trios[0][0].sample_id == "CHILD1"


def test_headered(tmp_path):
    ped = parse_ped(_write(tmp_path, HEADERED))
    assert len(ped.individuals) == 3
    assert len(ped.trios()) == 1
    assert ped.sex_of("B00ID5W") == 2
    assert ped.individuals["C0005IO"].phenotype == "-9"


def test_conflicting_duplicate_id_raises(tmp_path):
    content = HEADERLESS + "FAM9\tCHILD1\t0\t0\t1\t1\n"  # different family: conflict
    with pytest.raises(ValueError, match="duplicate"):
        parse_ped(_write(tmp_path, content))


def test_identical_duplicate_line_tolerated(tmp_path):
    content = HEADERLESS + "FAM1\tMOM1\t0\t0\t2\t1\n"  # exact repeat of MOM1
    ped = parse_ped(_write(tmp_path, content))
    assert len(ped.individuals) == 3


def test_duplicate_with_another_sex_flagged(tmp_path, caplog):
    content = HEADERLESS + "FAM1\tMOM1\t0\t0\t1\t1\n"  # MOM1 again, as a man
    ped = parse_ped(_write(tmp_path, content))
    assert ped.sex_of("MOM1") == 2                       # first line kept
    assert "another sex/phenotype" in caplog.text


def test_too_few_columns(tmp_path):
    with pytest.raises(ValueError, match="6"):
        parse_ped(_write(tmp_path, "A\tB\tC\n"))


def test_sex_written_as_float(tmp_path):
    from svengers.ped import parse_ped

    p = tmp_path / "ped.tsv"
    p.write_text("F1\tA\t0\t0\t1.0\t1\nF1\tB\t0\t0\t2.0\t1\nF1\tC\tA\tB\t-9\t2\n")
    ped = parse_ped(p)
    assert [ped.sex_of(s) for s in "ABC"] == [1, 2, 0]
