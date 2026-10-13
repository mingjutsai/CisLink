"""CisLink models: L2-regularized logistic regressions stored as plain JSON.

A model scores a table of harmonized features (one row per element-gene or variant-gene pair):

    score = sigmoid(intercept + coef . standardize(sign(x) * log1p(|x|)))

The JSON files hold the feature order, the standardization and the coefficients of the fitted
scikit-learn pipelines, so scoring needs only NumPy and gives the same values.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import resources
from typing import Iterable

import numpy as np
import pandas as pd

MODELS = ("e2g", "v2g")
# Hi-C contact is missing where a map has no record; the models were trained with these filled by 0.
FILL_ZERO = ("HiC_Contact", "HiC_megamap_oe")
# earlier column names still accepted (accessibility was called DNase, though ATAC-seq works as well)
LEGACY = {"Accessibility_Enh": "DNase_Enh", "Accessibility_Prom": "DNase_Prom"}


def slog(x: np.ndarray) -> np.ndarray:
    """The input transform of both models."""
    return np.sign(x) * np.log1p(np.abs(x))


@dataclass(frozen=True)
class CisLinkModel:
    name: str
    features: tuple[str, ...]
    scaler_mean: np.ndarray
    scaler_scale: np.ndarray
    coef: np.ndarray
    intercept: float
    thresholds: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @classmethod
    def from_json(cls, path_or_dict) -> "CisLinkModel":
        d = path_or_dict if isinstance(path_or_dict, dict) else json.load(open(path_or_dict))
        return cls(name=d["name"], features=tuple(d["features"]), scaler_mean=np.asarray(d["scaler_mean"], float),
                   scaler_scale=np.asarray(d["scaler_scale"], float), coef=np.asarray(d["coef"], float),
                   intercept=float(d["intercept"]), thresholds=d.get("thresholds", {}), meta=d)

    def matrix(self, table: pd.DataFrame) -> np.ndarray:
        """The model's input matrix from a feature table, with the checks applied before scoring."""
        table = table.rename(columns={old: new for new, old in LEGACY.items() if new not in table.columns and old in table.columns})
        missing = [f for f in self.features if f not in table.columns]
        if missing:
            raise KeyError(f"{self.name}: missing feature columns {missing}")
        x = table.loc[:, list(self.features)].copy()
        for f in FILL_ZERO:
            if f in x:
                x[f] = x[f].fillna(0.0)
        x = x.apply(pd.to_numeric, errors="raise")
        bad = [f for f in self.features if x[f].isna().any()]
        if bad:
            raise ValueError(f"{self.name}: missing values in {bad} (only Hi-C contact may be missing)")
        return x.to_numpy(dtype=float)

    def decision_function(self, table: pd.DataFrame) -> np.ndarray:
        z = (slog(self.matrix(table)) - self.scaler_mean) / self.scaler_scale
        return z @ self.coef + self.intercept

    def predict(self, table: pd.DataFrame) -> np.ndarray:
        """Scores in [0, 1]: probabilities on the scale of the K562 CRISPR training pairs."""
        return 1.0 / (1.0 + np.exp(-self.decision_function(table)))

    def contributions(self, table: pd.DataFrame) -> pd.DataFrame:
        """Per-feature contributions to the log-odds (coefficient x standardized input); rows sum to the
        log-odds minus the intercept."""
        z = (slog(self.matrix(table)) - self.scaler_mean) / self.scaler_scale
        return pd.DataFrame(z * self.coef, columns=list(self.features), index=table.index)

    def threshold(self, recall: int = 70) -> float:
        """Score at which held-out K562 CRISPR pairs reach the given recall (70 or 50)."""
        return float(self.thresholds[f"recall_{recall}"])


def load_model(name: str) -> CisLinkModel:
    """Load a packaged model: 'e2g' (CisLink-E2G) or 'v2g' (CisLink-V2G)."""
    key = name.lower().replace("cislink-", "")
    if key not in MODELS:
        raise ValueError(f"unknown model {name!r}; choose from {MODELS}")
    with resources.files("cislink.models").joinpath(f"{key}.json").open() as fh:
        return CisLinkModel.from_json(json.load(fh))


def score(table: pd.DataFrame, models: Iterable[str] = MODELS) -> pd.DataFrame:
    """Add a score column per model ('CisLink-E2G', 'CisLink-V2G') to a copy of the table."""
    out = table.copy()
    for m in models:
        model = load_model(m)
        out[model.name] = model.predict(table)
    return out
