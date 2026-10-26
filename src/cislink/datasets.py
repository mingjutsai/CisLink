"""Loaders for the CRISPR benchmark tables in ``data/crispr`` (repository checkout; not in the wheel)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .model import MODELS, load_model

CELLS = ("K562", "HCT116", "GM12878", "WTC11")   # WTC11 has no histone marks: CisLink-V2G only


def data_dir(path: str | Path | None = None) -> Path:
    if path is not None:
        return Path(path)
    for p in (Path.cwd(), *Path.cwd().parents):
        if (p / "data" / "crispr").is_dir():
            return p / "data" / "crispr"
    raise FileNotFoundError("data/crispr not found; run from the repository or pass data_dir")


def crispr_features(cell: str, path=None) -> pd.DataFrame:
    """Feature rows (one per element and transcript of the gene) of a held-out CRISPR set."""
    return pd.read_csv(data_dir(path) / f"{cell}_heldout.features.tsv.gz", sep="\t", low_memory=False)


def crispr_pairs(cell: str, path=None, score: bool = True) -> pd.DataFrame:
    """Element-gene pairs of a held-out CRISPR set with labels and comparator scores. With ``score``, adds
    CisLink-E2G and CisLink-V2G computed from the feature rows (maximum over a pair's transcripts)."""
    p = pd.read_csv(data_dir(path) / f"{cell}_heldout.pairs.tsv.gz", sep="\t")
    if score:
        f = crispr_features(cell, path)
        for m in MODELS:
            model = load_model(m)
            if model.missing(f.columns):   # e.g. CisLink-E2G in WTC11 (no histone marks)
                p[model.name] = float("nan")
                continue
            p[model.name] = p.ID.map(pd.Series(model.predict(f), index=f.ID).groupby(level=0).max()).fillna(0.0).values
    p["cluster"] = p.gene
    return p

