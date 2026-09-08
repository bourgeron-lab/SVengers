"""Tier-2 calibrated model scoring (SQ).

A trained model directory contains:
    model_del.joblib / model_dup.joblib   sklearn estimators with predict_proba
    manifest.json                         {model_version, schema_hash, features,
                                           calibration, training_data, thresholds}

Scoring refuses to run on a feature-schema mismatch. Until a trained model
ships, `vedette filter` runs Tier-1 only and leaves SQ absent.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path

import numpy as np

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

    @classmethod
    def load(cls, model_dir: str | Path) -> "FilterModel":
        import joblib

        model_dir = Path(model_dir)
        manifest = json.loads((model_dir / "manifest.json").read_text())
        if manifest["schema_hash"] != schema_hash():
            raise ValueError(
                f"Model {model_dir} was trained on feature schema "
                f"{manifest['schema_hash']}, this build has {schema_hash()}: "
                "retrain or use a matching vedette version"
            )
        return cls(
            joblib.load(model_dir / "model_del.joblib"),
            joblib.load(model_dir / "model_dup.joblib"),
            manifest,
        )

    def score(self, features: dict[str, float], svtype: str) -> float:
        """Phred-scaled calibrated P(true call)."""
        x = np.array(
            [[features[name] for name in FEATURES]], dtype=np.float64
        )
        model = self.model_dup if svtype == "DUP" else self.model_del
        p_true = float(model.predict_proba(x)[0, 1])
        p_true = min(max(p_true, 1e-10), 1 - 10 ** (-MAX_SQ / 10))
        return -10.0 * math.log10(1.0 - p_true)
