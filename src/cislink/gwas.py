"""Locus-level evaluation of variant-to-gene methods at GWAS loci (Supplementary Note 1).

A locus (credible set) has candidate genes with one score column per method and a boolean target column.
* Top-1 accuracy: the share of targets among the genes tied for first, averaged over loci. A method whose highest
  score at a locus is <= 0 makes no prediction there (top-1 = 0), except methods that always call (nearest gene).
* Reciprocal rank uses the best-ranked target (ties averaged); the chance-corrected version is (RR - E) / (1 - E),
  with E the expected reciprocal rank under random ordering.
* Locus-level precision-recall: each method calls its top gene, ordered by its top score; tied confidence blocks
  are drawn at their expected precision.
* 95% CIs: bootstrap over fine-mapping regions, so that credible sets of one region are resampled together.
"""
from __future__ import annotations

from math import comb

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def expected_rr(n: int, g: int) -> float:
    """Expected reciprocal rank of the best of ``g`` targets among ``n`` candidates in random order."""
    return sum((1 / k) * comb(n - k, g - 1) / comb(n, g) for k in range(1, n - g + 2))


def per_locus(t: pd.DataFrame, score: str, *, locus: str = "locus", target: str = "target", region: str = "region",
              always: bool = False) -> pd.DataFrame:
    rows = []
    for loc, s in t.groupby(locus, sort=False):
        x = s[score].fillna(-np.inf if always else 0).to_numpy(float); y = s[target].to_numpy(bool)
        mx = x.max(); top = x == mx
        t1 = 0.0 if (mx <= 0 and not always) else float(y[top].mean())
        best = pd.Series(x).rank(ascending=False, method="average").to_numpy()[y].min()
        e = expected_rr(len(s), int(y.sum()))
        rows.append({"locus": loc, "region": s[region].iloc[0], "top1": t1, "rr": 1 / best,
                     "cc": (1 / best - e) / (1 - e) if e < 1 else 0.0, "called": bool(mx > 0 or always)})
    return pd.DataFrame(rows)


def summarize(t: pd.DataFrame, methods: dict[str, bool], **kw) -> pd.DataFrame:
    """Top-1, MRR, chance-corrected MRR and share of loci called per method. ``methods``: {column: always_calls}."""
    out = []
    for m, always in methods.items():
        p = per_locus(t, m, always=always, **kw)
        out.append({"method": m, "loci": len(p), "top1": p.top1.mean(), "MRR": p.rr.mean(), "MRR_cc": p.cc.mean(), "called": p.called.mean()})
    return pd.DataFrame(out)


def region_bootstrap_difference(a: pd.DataFrame, b: pd.DataFrame, metric: str = "top1", n: int = 1000, seed: int = 0) -> tuple[float, float, float]:
    """Paired difference of a per-locus metric (outputs of ``per_locus`` on the same loci), bootstrap over regions."""
    A = a.set_index("locus"); Bm = b.set_index("locus").reindex(A.index)
    grp = {r: A.index[A.region == r] for r in A.region.unique()}; ids = list(grp); rng = np.random.default_rng(seed); bs = []
    for _ in range(n):
        sel = np.concatenate([grp[i] for i in rng.choice(ids, len(ids))])
        bs.append(A.loc[sel, metric].mean() - Bm.loc[sel, metric].mean())
    return A[metric].mean() - Bm[metric].mean(), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def locus_calls(t: pd.DataFrame, score: str, *, locus: str = "locus", target: str = "target") -> pd.DataFrame:
    """Per locus: confidence (top score; 0 = no call) and expected correctness of the top call (ties split)."""
    rows = []
    for _, s in t.groupby(locus, sort=False):
        x = s[score].fillna(0).to_numpy(float); mx = x.max()
        rows.append((mx if mx > 0 else 0.0, float(s[target].to_numpy(bool)[x == mx].mean()) if mx > 0 else 0.0))
    return pd.DataFrame(rows, columns=["conf", "correct"])


def locus_pr(calls: pd.DataFrame) -> dict:
    """Average precision and the precision-recall curve of locus calls."""
    L = len(calls); u = calls[calls.conf > 0]
    g = u.groupby("conf", sort=True).agg(n=("correct", "size"), k=("correct", "sum")).iloc[::-1]
    prec, rec = g.k.cumsum().to_numpy() / g.n.cumsum().to_numpy(), g.k.cumsum().to_numpy() / L
    return {"AP": float(np.sum(np.diff(np.concatenate([[0], rec])) * prec)), "precision": prec, "recall": rec,
            "max_recall": float(rec[-1]) if len(rec) else 0.0}


def pair_level(t: pd.DataFrame, score: str, *, nearest: str, locus: str = "locus", target: str = "target") -> dict:
    """Gene-locus pairs (L2G format): average precision and AUROC of ``score`` and of the naive nearest-gene classifier
    (1 for the gene nearest the lead variant, 0 otherwise)."""
    y = t[target].to_numpy(bool)
    naive = (t[nearest] == t.groupby(locus)[nearest].transform("max")).to_numpy(float)
    s = t[score].fillna(0).to_numpy(float)
    return {"AP": average_precision_score(y, s), "AUROC": roc_auc_score(y, s),
            "AP nearest": average_precision_score(y, naive), "AUROC nearest": roc_auc_score(y, naive)}
