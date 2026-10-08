"""Training of the calibrated DEL/DUP filter models.

Design (locked, see docs/ROADMAP.md):
- two separate HistGradientBoostingClassifier models (DEL, DUP);
- monotone constraints where the physics is unambiguous (e.g. more split
  reads can never make a call less credible);
- leave-chromosome-groups-out CV -> out-of-fold probabilities;
- isotonic calibration fitted on OOF probabilities (so SQ is a real,
  frequency-calibrated P(true));
- final base model refit on all data, shipped as CalibratedModel with a
  manifest pinning the feature-schema hash and threshold->FDR table.
"""

from __future__ import annotations

import csv
import json
import logging
import time
from pathlib import Path

import numpy as np

from .. import __version__
from ..filter.model import CalibratedModel
from ..filter.schema import FEATURES, schema_hash

logger = logging.getLogger(__name__)

# +1: P(true) non-decreasing in the feature; -1: non-increasing; 0: free
_MONOTONE = {
    "DEL": {"dhffc": -1, "fc_chrom": -1, "sr": +1, "pe": +1, "cig": +1, "ins_z": +1,
            "frac_consistent": +1},
    "DUP": {"dhffc": +1, "fc_chrom": +1, "sr": +1, "pe": +1, "frac_consistent": +1},
}

_HGB_PARAMS = dict(
    max_iter=400,
    learning_rate=0.06,
    max_leaf_nodes=31,
    min_samples_leaf=50,
    l2_regularization=1.0,
    random_state=17,
)


def load_labeled(paths: list[str | Path]) -> dict[str, dict[str, np.ndarray]]:
    """Load labeled feature TSVs -> per-svtype {X, y, contig}."""
    rows: dict[str, list] = {"DEL": [], "DUP": []}
    contigs: dict[str, list] = {"DEL": [], "DUP": []}
    labels: dict[str, list] = {"DEL": [], "DUP": []}
    for path in paths:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            missing = set(FEATURES + ["label", "svtype", "contig"]) - set(
                reader.fieldnames or ()
            )
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            for row in reader:
                svtype = row["svtype"]
                if svtype not in rows:
                    continue
                rows[svtype].append(
                    [float(row[f]) if row[f] not in ("", "nan") else np.nan for f in FEATURES]
                )
                contigs[svtype].append(row["contig"])
                labels[svtype].append(int(row["label"]))
    out = {}
    for svtype in ("DEL", "DUP"):
        out[svtype] = {
            "X": np.asarray(rows[svtype], dtype=np.float64),
            "y": np.asarray(labels[svtype], dtype=np.int64),
            "contig": np.asarray(contigs[svtype]),
        }
        logger.info(
            "%s: %d calls (%d TP / %d FP)",
            svtype, len(labels[svtype]),
            int(out[svtype]["y"].sum()),
            int((out[svtype]["y"] == 0).sum()),
        )
    return out


# Features that are legitimately constant within one SV type's training set:
# `is_dup` because DEL and DUP get their own model, and the DEL-specific
# evidence channels, which are structurally zero on duplications.
_EXPECTED_CONSTANT: dict[str, set[str]] = {
    "DEL": {"is_dup"},
    "DUP": {"is_dup", "cig", "ins_z", "ins_pairs"},
}


def constant_features(X: np.ndarray, svtype: str, skip: set[str] | frozenset = frozenset()) -> dict[str, str]:
    """Unexpectedly constant (or all-missing) feature columns -> why.

    A feature with no variance carries no information: the trees never split
    on it and it silently contributes nothing. `lcr_ovl` sat dead like this
    through every model up to v5 — the GHFC LCR track is Ensembl-named while
    calls are UCSC, so the overlap was always 0 (see docs/ROADMAP.md) — and
    nothing in the pipeline complained. This guard makes that fail loudly.
    """
    expected = _EXPECTED_CONSTANT.get(svtype, {"is_dup"}) | set(skip)
    dead: dict[str, str] = {}
    for i, name in enumerate(FEATURES):
        if name in expected:
            continue
        col = X[:, i]
        finite = col[np.isfinite(col)]
        if len(finite) == 0:
            dead[name] = "aucune valeur finie"
        elif len(np.unique(finite)) == 1:
            dead[name] = f"constante = {finite[0]:g}"
    return dead


def check_feature_variance(
    X: np.ndarray, svtype: str, allow_constant: bool = False,
    skip: set[str] | frozenset = frozenset(),
) -> dict[str, str]:
    """Raise on unexpectedly constant features unless explicitly allowed."""
    dead = constant_features(X, svtype, skip)
    if not dead:
        return dead
    detail = ", ".join(f"{k} ({v})" for k, v in sorted(dead.items()))
    if allow_constant:
        logger.warning("%s: feature(s) sans variance, entraînement forcé: %s", svtype, detail)
        return dead
    raise ValueError(
        f"{svtype}: {len(dead)} feature(s) sans variance dans le jeu d'entraînement: "
        f"{detail}. Une feature constante n'apporte rien au modèle et trahit "
        f"en général une extraction cassée (piste: convention de nommage des "
        f"contigs du track, chemin de ressource absent). Corrigez l'extraction, "
        f"ou passez --allow-constant-features si c'est intentionnel."
    )


def chromosome_folds(contigs: np.ndarray, n_folds: int = 5) -> np.ndarray:
    """Assign whole chromosomes to folds, balancing call counts (greedy)."""
    uniq, counts = np.unique(contigs, return_counts=True)
    order = np.argsort(-counts)
    fold_of: dict[str, int] = {}
    load = np.zeros(n_folds)
    for idx in order:
        f = int(np.argmin(load))
        fold_of[uniq[idx]] = f
        load[f] += counts[idx]
    return np.asarray([fold_of[c] for c in contigs])


def _fit_base(X: np.ndarray, y: np.ndarray, svtype: str):
    from sklearn.ensemble import HistGradientBoostingClassifier

    cst = [_MONOTONE[svtype].get(name, 0) for name in FEATURES]
    model = HistGradientBoostingClassifier(monotonic_cst=cst, **_HGB_PARAMS)
    model.fit(X, y)
    return model


def train_one(
    X: np.ndarray, y: np.ndarray, contigs: np.ndarray, svtype: str, n_folds: int = 5
) -> tuple[CalibratedModel, dict]:
    from sklearn.isotonic import IsotonicRegression
    from sklearn.metrics import average_precision_score, roc_auc_score

    if len(np.unique(y)) < 2:
        raise ValueError(f"{svtype}: need both classes to train (got only {set(y)})")
    folds = chromosome_folds(contigs, n_folds)
    oof = np.full(len(y), np.nan)
    for f in range(n_folds):
        test = folds == f
        if not test.any() or len(np.unique(y[~test])) < 2:
            continue
        model = _fit_base(X[~test], y[~test], svtype)
        oof[test] = model.predict_proba(X[test])[:, 1]
    valid = ~np.isnan(oof)

    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(oof[valid], y[valid])
    calibrated_oof = iso.predict(oof[valid])

    base = _fit_base(X, y, svtype)
    model = CalibratedModel(base, iso)

    # threshold -> achieved FDR / recall table on calibrated OOF probabilities
    yv = y[valid]
    table = []
    for sq in (3, 7, 10, 13, 17, 20, 30):
        p_min = 1 - 10 ** (-sq / 10)
        kept = calibrated_oof >= p_min
        if kept.sum() == 0:
            table.append({"sq": sq, "kept": 0, "fdr": None, "recall": 0.0})
            continue
        fdr = float((yv[kept] == 0).mean())
        recall = float(yv[kept].sum() / max(yv.sum(), 1))
        table.append({"sq": sq, "kept": int(kept.sum()), "fdr": round(fdr, 4),
                      "recall": round(recall, 4)})

    metrics = {
        "n": int(len(y)),
        "n_tp": int(y.sum()),
        "n_labeled_oof": int(valid.sum()),
        "pr_auc_oof": round(float(average_precision_score(yv, oof[valid])), 4),
        "roc_auc_oof": round(float(roc_auc_score(yv, oof[valid])), 4),
        "threshold_table": table,
    }
    logger.info("%s: PR-AUC(OOF)=%.4f ROC-AUC=%.4f",
                svtype, metrics["pr_auc_oof"], metrics["roc_auc_oof"])
    for row in table:
        logger.info("%s  SQ>=%d: kept=%d FDR=%s recall=%s",
                    svtype, row["sq"], row["kept"], row["fdr"], row["recall"])
    return model, metrics


def train_models(
    labeled_tsvs: list[str | Path],
    out_dir: str | Path,
    model_version: str = "v1",
    training_data_desc: str = "",
    allow_constant_features: bool = False,
    neutral_features: dict[str, float] | None = None,
) -> dict:
    """`neutral_features` ({name: value}) are pinned to a constant here and,
    through the manifest, at scoring (see filter/model.py)."""
    import joblib

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = load_labeled(labeled_tsvs)
    neutral = dict(neutral_features or {})
    unknown = set(neutral) - set(FEATURES)
    if unknown:
        raise ValueError(f"neutral features absent from the schema: {sorted(unknown)}")
    for d in data.values():
        for name, value in neutral.items():
            d["X"][:, FEATURES.index(name)] = value
    manifest = {
        "model_version": model_version,
        "vedette_version": __version__,
        "schema_hash": schema_hash(),
        "features": FEATURES,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "training_data": training_data_desc,
        "metrics": {},
    }
    if neutral:
        manifest["neutral_features"] = neutral
    # fail before any expensive fitting if a feature carries no signal
    for svtype in ("DEL", "DUP"):
        dead = check_feature_variance(
            data[svtype]["X"], svtype, allow_constant=allow_constant_features,
            skip=set(neutral),
        )
        if dead:
            manifest.setdefault("constant_features", {})[svtype] = sorted(dead)
    for svtype, fname in (("DEL", "model_del.joblib"), ("DUP", "model_dup.joblib")):
        d = data[svtype]
        model, metrics = train_one(d["X"], d["y"], d["contig"], svtype)
        joblib.dump(model, out_dir / fname)
        manifest["metrics"][svtype] = metrics
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    logger.info("models written to %s", out_dir)
    return manifest
