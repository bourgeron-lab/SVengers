"""Tier-2 calibrated model scoring (SQ).

A trained model directory contains:
    model_del.joblib / model_dup.joblib   sklearn estimators with predict_proba
    manifest.json                         {model_version, schema_hash, features,
                                           calibration, training_data, thresholds,
                                           neutral_features}

`neutral_features` ({name: value}) pins features to a constant at training
and scoring alike. `wave_score` is pinned in v5n/v8: v5 learned a step on
it between 0.10 and 0.13 — the same calls of a PMS genome went from 4,794
to 29,444 at SQ >= 5 with wave 0.10 vs 0.13 — so every genome past the step
was flooded with PASS calls (a third of PMS). dataset_paper escaped only
because its stores predate the feature (-1).

Scoring refuses to run on a feature-schema mismatch. Until a trained model
ships, `svengers filter` runs Tier-1 only and leaves SQ absent.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path

import numpy as np

from ..legacy import caller_version, legacy_modules
from .schema import FEATURES, schema_hash

logger = logging.getLogger(__name__)

MAX_SQ = 99.0


class CalibratedModel:
    """A fitted classifier + isotonic calibrator, exposing predict_proba.

    Defined here (not in train/) so joblib artifacts stay loadable in
    runtime-only installs.
    """

    def __init__(self, base, calibrator):
        self.base = base
        self.calibrator = calibrator

    def predict_proba(self, X) -> np.ndarray:
        raw = self.base.predict_proba(X)[:, 1]
        cal = np.clip(self.calibrator.predict(raw), 0.0, 1.0)
        return np.stack([1.0 - cal, cal], axis=1)


class FilterModel:
    def __init__(self, model_del, model_dup, manifest: dict):
        self.model_del = model_del
        self.model_dup = model_dup
        self.manifest = manifest
        self.version = manifest.get("model_version", "unknown")
        # caller version that trained it (vedette_version in pre-2.0 manifests)
        self.caller_version = caller_version(manifest)
        self.features = manifest.get("features") or list(FEATURES)
        self.neutral = {
            self.features.index(name): float(value)
            for name, value in (manifest.get("neutral_features") or {}).items()
            if name in self.features
        }

    def _neutralize(self, x: np.ndarray) -> np.ndarray:
        for i, value in self.neutral.items():
            x[:, i] = value
        return x

    @classmethod
    def load(cls, model_dir: str | Path) -> "FilterModel":
        import joblib

        model_dir = Path(model_dir)
        manifest = json.loads((model_dir / "manifest.json").read_text())
        if manifest["schema_hash"] != schema_hash():
            # append-only contract: a model trained on a PREFIX of the current
            # schema still scores (with its own feature list); anything else
            # is a real mismatch
            trained = manifest.get("features")
            if not trained or trained != FEATURES[: len(trained)]:
                raise ValueError(
                    f"Model {model_dir} was trained on feature schema "
                    f"{manifest['schema_hash']}, this build has {schema_hash()}: "
                    "retrain or use a matching svengers version"
                )
        with legacy_modules():  # models pickled before 2.0 name vedette.filter.model
            return cls(
                joblib.load(model_dir / "model_del.joblib"),
                joblib.load(model_dir / "model_dup.joblib"),
                manifest,
            )

    def score(self, features: dict[str, float], svtype: str) -> float:
        """Phred-scaled calibrated P(true call)."""
        x = self._neutralize(np.array(
            [[features[name] for name in self.features]], dtype=np.float64
        ))
        model = self.model_dup if svtype == "DUP" else self.model_del
        p_true = float(model.predict_proba(x)[0, 1])
        p_true = min(max(p_true, 1e-10), 1 - 10 ** (-MAX_SQ / 10))
        return -10.0 * math.log10(1.0 - p_true)

    def score_batch(
        self, feature_rows: list[dict[str, float]], svtypes: list[str]
    ) -> np.ndarray:
        """Vectorized `score` over many calls (one predict_proba per SV type).

        Per-call predict_proba dominates the filter runtime with deep forests
        (~5 h/genome with v5); one batched pass is ~100x faster and returns
        bit-identical SQ values.
        """
        x = self._neutralize(np.array(
            [[row[name] for name in self.features] for row in feature_rows],
            dtype=np.float64,
        ))
        types = np.array(svtypes)
        sq = np.zeros(len(feature_rows), dtype=np.float64)
        cap = 1 - 10 ** (-MAX_SQ / 10)
        for svtype, model in (("DUP", self.model_dup), ("DEL", self.model_del)):
            mask = (types == "DUP") if svtype == "DUP" else (types != "DUP")
            if mask.any():
                p = model.predict_proba(x[mask])[:, 1]
                p = np.clip(p, 1e-10, cap)
                sq[mask] = -10.0 * np.log10(1.0 - p)
        return sq
