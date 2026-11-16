"""SVengers was Vedette up to v1.16: what Vedette wrote stays readable."""

from __future__ import annotations

import json
import re
import shutil
import tomllib
from pathlib import Path

import pysam
import pytest

from svengers.export import family_export_path, read_family_meta
from svengers.filter.model import FilterModel
from svengers.legacy import caller_version, existing_raw_calls, raw_calls_vcf
from svengers.vcfout import add_coords_line, span, spec_coords

ROOT = Path(__file__).resolve().parents[1]


def _header(*lines: str) -> pysam.VariantHeader:
    h = pysam.VariantHeader()
    for line in lines:
        h.add_line(line)
    return h


def test_spec_coords_both_names():
    # padded coordinates, marked with either name
    assert spec_coords(_header("##source=vedette-1.16.0", "##vedette_coords=padded"))
    assert spec_coords(_header("##source=svengers-2.0.0", "##svengers_coords=padded"))
    # a padded Vedette file run through another tool that changed the source line
    assert spec_coords(_header("##source=bcftools", "##vedette_coords=padded"))
    # our own source without the coordinates line: pre-1.4 layout, either name
    assert not spec_coords(_header("##source=vedette-1.3.3"))
    assert not spec_coords(_header("##source=vedette-merge-1.3.3"))
    assert not spec_coords(_header("##source=svengers-1.3.3"))
    # any other caller follows the spec
    assert spec_coords(_header("##source=Manta"))


def test_add_coords_line_keeps_a_legacy_key():
    h = _header("##source=vedette-1.16.0", "##vedette_coords=padded")
    add_coords_line(h)
    assert sum(r.key.endswith("_coords") for r in h.records) == 1
    h = _header("##source=svengers-2.0.0")
    add_coords_line(h)
    assert [r.value for r in h.records if r.key == "svengers_coords"] == ["padded"]


def test_padded_vedette_vcf_read_as_padded(tmp_path):
    h = _header("##source=vedette-1.16.0", "##vedette_coords=padded")
    h.contigs.add("chr1", length=1_000_000)
    h.add_line('##ALT=<ID=DEL,Description="Deletion">')
    h.info.add("SVTYPE", "1", "String", "x")
    h.info.add("SVLEN", "1", "Integer", "x")
    path = tmp_path / "S1.vedette.vcf"
    with pysam.VariantFile(str(path), "w", header=h) as o:
        rec = o.new_record(contig="chr1", start=999, stop=1_100, alleles=("N", "<DEL>"))
        rec.info["SVTYPE"], rec.info["SVLEN"] = "DEL", -100
        o.write(rec)
    with pysam.VariantFile(str(path)) as vf:
        assert spec_coords(vf.header)
        (rec,) = vf
        assert span(rec, spec_coords(vf.header)) == (1_000, 1_100)


def test_raw_calls_lookup_falls_back_to_vedette_name(tmp_path):
    d = tmp_path / "SMP-01"
    d.mkdir()
    # neither: the new name (for the error message)
    assert raw_calls_vcf(d, "SMP-01") == d / "SMP-01.svengers.vcf.gz"
    old = d / "SMP-01.vedette.vcf.gz"
    old.touch()
    assert raw_calls_vcf(d, "SMP-01") == old
    assert existing_raw_calls(d / "SMP-01.svengers.vcf.gz") == old
    # the new name wins when both exist
    new = d / "SMP-01.svengers.vcf.gz"
    new.touch()
    assert raw_calls_vcf(d, "SMP-01") == new
    assert existing_raw_calls(new) == new
    # any other path is left alone
    other = d / "SMP-01.filtered.vcf.gz"
    assert existing_raw_calls(other) == other


def test_caller_version_both_keys():
    assert caller_version({"svengers_version": "2.0.0"}) == "2.0.0"
    assert caller_version({"vedette_version": "1.16.0"}) == "1.16.0"
    assert caller_version({"svengers_version": "2.0.0", "vedette_version": "1.16.0"}) == "2.0.0"
    assert caller_version({}) is None


def test_read_family_meta_of_a_vedette_export(tmp_path):
    info = {"schema_version": 1, "vedette_version": "1.16.0", "family_id": "FAM-01",
            "members": ["SMP-01", "SMP-02", "SMP-03"]}
    (tmp_path / "FAM-01.vedette.json").write_text(json.dumps(info))
    assert family_export_path(tmp_path, "FAM-01", "json") == tmp_path / "FAM-01.vedette.json"
    assert family_export_path(tmp_path, "FAM-01", "tsv") == tmp_path / "FAM-01.svengers.tsv"
    meta = read_family_meta(tmp_path, "FAM-01")
    assert meta["svengers_version"] == meta["vedette_version"] == "1.16.0"
    assert meta["members"] == info["members"]


def test_read_family_meta_from_a_vedette_parquet(tmp_path):
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")
    info = {"schema_version": 1, "vedette_version": "1.15.1", "family_id": "FAM-01"}
    table = pa.table({"row_id": ["x"]}).replace_schema_metadata(
        {"vedette": json.dumps(info), "vedette_version": "1.15.1"})
    pq.write_table(table, tmp_path / "FAM-01.vedette.parquet")
    assert read_family_meta(tmp_path, "FAM-01")["svengers_version"] == "1.15.1"


def test_model_manifest_with_vedette_version(tmp_path):
    # a tracked model trained before the rename loads, version read back
    tracked = ROOT / "models" / "v12b"
    model = FilterModel.load(tracked)
    manifest = json.loads((tracked / "manifest.json").read_text())
    assert "vedette_version" in manifest and "svengers_version" not in manifest
    assert model.caller_version == manifest["vedette_version"]
    # same manifest with the new key
    copy = tmp_path / "m"
    shutil.copytree(tracked, copy)
    manifest["svengers_version"] = manifest.pop("vedette_version")
    (copy / "manifest.json").write_text(json.dumps(manifest))
    assert FilterModel.load(copy).caller_version == model.caller_version


def test_both_console_scripts_declared():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["name"] == "svengers"
    assert project["scripts"] == {"svengers": "svengers.cli:app", "vedette": "svengers.cli:app"}


def test_nextflow_accepts_old_param_names_and_call_names():
    nf = ROOT / "workflows" / "nextflow"
    config = (nf / "nextflow.config").read_text()
    assert re.search(r"svengers_bin\s*=\s*params\.vedette_bin\s*\?:", config)
    assert re.search(r"svengers_config\s*=\s*params\.vedette_config\s*\?:", config)
    main = (nf / "main.nf").read_text()
    assert '"${dir}/${sample}.vedette.vcf.gz"' in main
