"""Command line: ``cislink score`` and ``cislink info``."""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from . import __version__
from .model import MODELS, load_model


def _score(a) -> int:
    t = pd.read_csv(a.input, sep=a.sep, low_memory=False)
    for m in a.model:
        model = load_model(m)
        t[model.name] = model.predict(t)
        if a.threshold:
            t[f"{model.name} >= {a.threshold}% recall threshold"] = t[model.name] >= model.threshold(a.threshold)
    keep = [c for c in a.keep.split(",") if c] if a.keep else None
    if keep:
        t = t[keep + [c for c in t.columns if c.startswith("CisLink-")]]
    t.to_csv(a.output or sys.stdout, sep="\t", index=False, float_format="%.6g")
    return 0


def _info(a) -> int:
    for m in MODELS:
        model = load_model(m)
        print(f"{model.name}: {len(model.features)} features; thresholds {model.thresholds.get('recall_70')} (70% recall), "
              f"{model.thresholds.get('recall_50')} (50% recall)")
        print("  " + ", ".join(model.features))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="cislink", description="Score element-gene and variant-gene pairs with CisLink.")
    p.add_argument("--version", action="version", version=f"cislink {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score", help="score a harmonized feature table")
    s.add_argument("input", help="feature table (TSV by default) with the model's feature columns")
    s.add_argument("-m", "--model", nargs="+", default=list(MODELS), choices=MODELS)
    s.add_argument("-o", "--output", help="output TSV (default stdout)")
    s.add_argument("--sep", default="\t")
    s.add_argument("--keep", help="comma-separated columns to keep with the scores (default: all)")
    s.add_argument("--threshold", type=int, choices=(50, 70), help="add a column flagging scores above the recall threshold")
    s.set_defaults(func=_score)
    i = sub.add_parser("info", help="list the models' features and thresholds")
    i.set_defaults(func=_info)
    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
