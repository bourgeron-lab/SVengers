import json

import numpy as np
import pytest

from vedette.filter.model import FilterModel
from vedette.filter.schema import FEATURES, schema_hash
from vedette.train.train import chromosome_folds, train_models

N_PER_TYPE = 3000


def _synthetic_labeled_tsv(tmp_path, seed=0):
    """Feature TSV where dhffc/sr separate TP from FP, per svtype physics."""
    rng = np.random.default_rng(seed)
    path = tmp_path / "labeled.tsv"
    header = ["call_id", "contig", "start", "end", "svtype", *FEATURES, "label"]
    contigs = [f"chr{i}" for i in range(1, 23)]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i in range(2 * N_PER_TYPE):
            svtype = "DEL" if i < N_PER_TYPE else "DUP"
            label = int(rng.random() < 0.4)
            feats = {name: 0.0 for name in FEATURES}
            if svtype == "DEL":
                feats["dhffc"] = rng.normal(0.55 if label else 0.95, 0.1)
            else:
                feats["dhffc"] = rng.normal(1.5 if label else 1.02, 0.12)
            feats["is_dup"] = float(svtype == "DUP")
            feats["sr"] = float(rng.poisson(8 if label else 0.4))
            feats["pe"] = float(rng.poisson(6 if label else 0.5))
            feats["size_log10"] = rng.uniform(2, 6)
            feats["baf_dev"] = -1.0
            if rng.random() < 0.05:
                feats["depth_cv"] = float("nan")  # NaN robustness
            row = [
                f"call_{i}", contigs[i % len(contigs)], str(1000 * i),
                str(1000 * i + 500), svtype,
                *(f"{feats[name]:.6g}" for name in FEATURES), str(label),
            ]
            fh.write("\t".join(row) + "\n")
    return path


def test_chromosome_folds_are_whole_chromosome_and_balanced():
    contigs = np.repeat([f"chr{i}" for i in range(1, 11)], 100)
    folds = chromosome_folds(contigs, n_folds=5)
    for c in np.unique(contigs):
        assert len(np.unique(folds[contigs == c])) == 1
    _, counts = np.unique(folds, return_counts=True)
    assert counts.min() >= 100 and counts.max() <= 300


def test_train_and_score_roundtrip(tmp_path):
    tsv = _synthetic_labeled_tsv(tmp_path)
    out = tmp_path / "models"
    manifest = train_models([tsv], out, model_version="test")

    for svtype in ("DEL", "DUP"):
        assert manifest["metrics"][svtype]["pr_auc_oof"] > 0.9
    saved = json.loads((out / "manifest.json").read_text())
    assert saved["schema_hash"] == schema_hash()

    model = FilterModel.load(out)
    good_del = {name: 0.0 for name in FEATURES}
    good_del.update(dhffc=0.5, sr=10.0, pe=8.0, size_log10=4.0, baf_dev=-1.0)
    bad_del = dict(good_del, dhffc=1.0, sr=0.0, pe=0.0)
    assert model.score(good_del, "DEL") > model.score(bad_del, "DEL")
    assert model.score(good_del, "DEL") > 10  # calibrated P(true) > 0.9

    good_dup = dict(good_del, dhffc=1.55, is_dup=1.0)
    bad_dup = dict(good_dup, dhffc=1.0, sr=0.0, pe=0.0)
    assert model.score(good_dup, "DUP") > model.score(bad_dup, "DUP")


def test_schema_mismatch_refused(tmp_path):
    # append-only contract: a manifest whose features are NOT a prefix of the
    # current schema (reordered here) is refused even if it self-declares a hash
    tsv = _synthetic_labeled_tsv(tmp_path)
    out = tmp_path / "models"
    train_models([tsv], out, model_version="test")
    manifest = json.loads((out / "manifest.json").read_text())
    manifest["schema_hash"] = "0" * 16
    manifest["features"] = manifest["features"][::-1]
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="schema"):
        FilterModel.load(out)


def test_prefix_schema_model_still_loads_and_scores(tmp_path):
    # a model trained before a feature was appended (its features are a strict
    # prefix of the current schema) keeps loading and scores on its own list
    tsv = _synthetic_labeled_tsv(tmp_path)
    out = tmp_path / "models"
    train_models([tsv], out, model_version="test")
    manifest = json.loads((out / "manifest.json").read_text())
    manifest["schema_hash"] = "0" * 16  # simulate: hashed before the append
    (out / "manifest.json").write_text(json.dumps(manifest))
    model = FilterModel.load(out)
    assert model.features == FEATURES  # trained on today's full list
    feats = {name: 0.0 for name in FEATURES}
    assert model.score(feats, "DEL") >= 0.0
