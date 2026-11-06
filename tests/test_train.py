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
            # every feature gets noise: a flat-zero fixture would trip the
            # zero-variance guard, which is exactly what it is there for
            feats = {name: float(rng.random()) for name in FEATURES}
            if svtype == "DEL":
                feats["dhffc"] = rng.normal(0.55 if label else 0.95, 0.1)
            else:
                feats["dhffc"] = rng.normal(1.5 if label else 1.02, 0.12)
            feats["is_dup"] = float(svtype == "DUP")
            feats["sr"] = float(rng.poisson(8 if label else 0.4))
            feats["pe"] = float(rng.poisson(6 if label else 0.5))
            feats["size_log10"] = rng.uniform(2, 6)
            # -1 is the "no het SNV here" sentinel; real cohorts mix both
            feats["baf_dev"] = -1.0 if rng.random() < 0.3 else float(rng.random() * 0.5)
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


def test_constant_feature_guard():
    """A feature with no variance must fail the training loudly: lcr_ovl sat
    dead (Ensembl-vs-UCSC contig naming) through every model up to v5 and
    nothing complained."""
    import numpy as np
    import pytest

    from vedette.filter.schema import FEATURES
    from vedette.train.train import check_feature_variance, constant_features

    rng = np.random.default_rng(0)
    X = rng.random((200, len(FEATURES)))
    i_is_dup = FEATURES.index("is_dup")
    i_lcr = FEATURES.index("lcr_ovl")
    i_cig = FEATURES.index("cig")

    # is_dup is constant by design (one model per SV type) -> not reported
    X[:, i_is_dup] = 0.0
    assert constant_features(X, "DEL") == {}
    assert check_feature_variance(X, "DEL") == {}

    # a genuinely dead feature is reported and raises
    X[:, i_lcr] = 0.0
    assert set(constant_features(X, "DEL")) == {"lcr_ovl"}
    with pytest.raises(ValueError, match="lcr_ovl"):
        check_feature_variance(X, "DEL")
    # ... unless explicitly allowed
    assert set(check_feature_variance(X, "DEL", allow_constant=True)) == {"lcr_ovl"}

    # DEL-only evidence channels are exempt on DUP but not on DEL
    X[:, i_cig] = 0.0
    X[:, i_is_dup] = 1.0
    assert set(constant_features(X, "DUP")) == {"lcr_ovl"}
    assert set(constant_features(X, "DEL")) == {"lcr_ovl", "cig"}

    # an all-NaN column counts as dead too
    X[:, i_lcr] = rng.random(200)
    X[:, i_cig] = np.nan
    assert constant_features(X, "DEL")["cig"] == "aucune valeur finie"


def test_neutral_feature_is_pinned_at_training_and_scoring(tmp_path):
    """v5 learned a step on the per-genome wave_score (same calls of a PMS
    genome: 4,794 vs 29,444 PASS-level calls at wave 0.10 vs 0.13). A neutral
    feature is constant in training and overridden at scoring, so a genome's
    wave can no longer move any call's SQ."""
    tsv = _synthetic_labeled_tsv(tmp_path)
    out = tmp_path / "models"
    manifest = train_models([tsv], out, model_version="test",
                            neutral_features={"wave_score": -1.0})
    assert manifest["neutral_features"] == {"wave_score": -1.0}
    model = FilterModel.load(out)
    call = {name: 0.3 for name in FEATURES}
    call.update(dhffc=0.6, sr=3.0, pe=2.0)
    scores = {w: model.score(dict(call, wave_score=w), "DEL") for w in (-1.0, 0.1, 0.16)}
    assert len(set(scores.values())) == 1
    batch = model.score_batch([dict(call, wave_score=w) for w in (-1.0, 0.16)], ["DEL", "DEL"])
    assert batch[0] == batch[1] == scores[-1.0]
    with pytest.raises(ValueError, match="neutral"):
        train_models([tsv], tmp_path / "bad", neutral_features={"nope": 0.0})


def test_weight_column_lets_a_few_curated_labels_count(tmp_path):
    """A few curated (IGV) labels against thousands of others: the optional
    ``weight`` column makes them count. Here 30 DEL calls that look like the
    ~1,800 false ones (flat depth, no junction support) are curated true: at
    weight 1 they are noise, at weight 100 they outweigh the region."""
    base = _synthetic_labeled_tsv(tmp_path)
    lines = base.read_text().splitlines()
    header = lines[0].split("\t")

    def write(path, weight):
        rng = np.random.default_rng(1)
        with open(path, "w") as fh:
            fh.write("\t".join([*header, "weight"]) + "\n")
            for line in lines[1:]:
                fh.write(line + "\t1\n")
            for i in range(30):
                feats = {name: float(rng.random()) for name in FEATURES}
                feats.update(dhffc=rng.normal(0.95, 0.1), sr=float(rng.poisson(0.4)),
                             pe=float(rng.poisson(0.5)), is_dup=0.0, size_log10=rng.uniform(2, 6))
                row = [f"igv_{i}", f"chr{1 + i % 22}", str(10**8 + 1000 * i), str(10**8 + 1000 * i + 500),
                       "DEL", *(f"{feats[name]:.6g}" for name in FEATURES), "1", str(weight)]
                fh.write("\t".join(row) + "\n")
        return path

    call = {name: 0.5 for name in FEATURES}
    call.update(dhffc=0.95, sr=0.0, pe=0.0, size_log10=3.5)
    scores = {}
    for weight in (1, 100):
        out = tmp_path / f"m{weight}"
        manifest = train_models([write(tmp_path / f"w{weight}.tsv", weight)], out)
        scores[weight] = FilterModel.load(out).score(call, "DEL")
        assert manifest["metrics"]["DEL"]["n"] == N_PER_TYPE + 30  # metrics count rows
    assert "weighted_rows" not in json.loads((tmp_path / "m1" / "manifest.json").read_text())["metrics"]["DEL"]
    assert manifest["metrics"]["DEL"]["weighted_rows"] == 30
    assert scores[1] < 1 and scores[100] > scores[1] + 2
