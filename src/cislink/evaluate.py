"""Evaluation used in the manuscript: AUPRC with gene-clustered bootstrap CIs, paired differences,
pooling across CRISPR sets, calibration, recall-based thresholds and grouped permutation importance."""
from __future__ import annotations

from typing import Callable, Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score

B, SEED = 1000, 0
DISTANCE_BINS = ((0, 1e4, "<10 kb"), (1e4, 5e4, "10-50 kb"), (5e4, 2e5, "50-200 kb"), (2e5, np.inf, ">200 kb"))
GROUPS = {
    "accessibility": ["Accessibility_Enh", "Accessibility_Prom", "DxA"],   # DxA = log10 distance x accessibility
    "histone marks": [f"{m}_{p}" for m in ("H4K20me1", "H3K9me3", "H3K9ac", "H3K79me2", "H3K4me3", "H3K4me2", "H3K4me1",
                                           "H3K36me3", "H3K27me3", "H3K27ac", "H2A.Z") for p in ("Enh", "Prom")],
    "CTCF": ["CTCF_Enh", "CTCF_Prom"],
    "expression": ["Tx_expression", "Gene_expression", "ubiquitous_gene"],
    "contact": ["HiC_megamap_oe", "HiC_Contact", "powerlaw_contact"],
    "distance and gene context": ["distance", "log10_distance", "n_genes_between", "gene_rank", "is_nearest_gene"],
    "element": ["Enh_length"],
    "cross-epigenome correlation": ["F1_dnase", "F1_h3k27ac"],
}


def auprc(y, s) -> float:
    return float(average_precision_score(np.asarray(y), np.nan_to_num(np.asarray(s, float))))


def cluster_bootstrap(clusters, b: int = B, seed: int = SEED, require=None) -> list[np.ndarray]:
    """Index sets resampling whole clusters (genes) with replacement. ``require``: a 0/1 label array; replicates
    with no positive are dropped."""
    cl = np.asarray(clusters)
    groups = pd.Series(np.arange(len(cl))).groupby(cl).apply(np.array).to_numpy()
    rng = np.random.default_rng(seed); out = []
    for _ in range(b):
        ii = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        if require is None or np.asarray(require)[ii].sum() > 0:
            out.append(ii)
    return out


def auprc_ci(y, s, clusters, b: int = B, seed: int = SEED) -> tuple[float, float, float]:
    y, s = np.asarray(y), np.asarray(s, float)
    bs = [auprc(y[i], s[i]) for i in cluster_bootstrap(clusters, b, seed, y)]
    return auprc(y, s), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def paired_difference(y, s1, s2, clusters, b: int = B, seed: int = SEED) -> tuple[float, float, float]:
    """AUPRC(s1) - AUPRC(s2) on the same pairs, with a gene-clustered bootstrap 95% CI."""
    y, s1, s2 = np.asarray(y), np.asarray(s1, float), np.asarray(s2, float)
    bs = [auprc(y[i], s1[i]) - auprc(y[i], s2[i]) for i in cluster_bootstrap(clusters, b, seed, y)]
    return auprc(y, s1) - auprc(y, s2), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def pooled_difference(sets: dict[str, pd.DataFrame], a: str, b_: str, b: int = B, seed: int = SEED) -> dict:
    """Positive-weighted mean of per-set AUPRC differences (columns ``a`` - ``b_``), bootstrapping genes within
    each set. Each set needs columns y, cluster, a and b_."""
    w = {k: int(t.y.sum()) for k, t in sets.items()}
    per = {k: auprc(t.y, t[a]) - auprc(t.y, t[b_]) for k, t in sets.items()}
    reps = {k: cluster_bootstrap(t.cluster.values, b, seed, t.y.values) for k, t in sets.items()}
    n = min(len(r) for r in reps.values()); vals = []
    for r in range(n):
        tot = 0.0
        for k, t in sets.items():
            i = reps[k][r]; tot += w[k] * (auprc(t.y.values[i], t[a].values[i]) - auprc(t.y.values[i], t[b_].values[i]))
        vals.append(tot / sum(w.values()))
    return {"difference": sum(per[k] * w[k] for k in per) / sum(w.values()), "lo": float(np.percentile(vals, 2.5)),
            "hi": float(np.percentile(vals, 97.5)), "per_set": per, "positives": w}


def distance_bins(distance) -> np.ndarray:
    d = np.abs(np.asarray(distance, float)); lab = np.empty(len(d), dtype=object)
    for lo, hi, name in DISTANCE_BINS:
        lab[(d >= lo) & (d < hi)] = name
    return lab


def calibration(y, s) -> dict:
    """Calibration slope and intercept: logistic regression of the label on logit(score)."""
    p = np.clip(np.asarray(s, float), 1e-4, 1 - 1e-4); z = np.log(p / (1 - p)).reshape(-1, 1)
    f = LogisticRegression(C=1e10, max_iter=1000).fit(z, np.asarray(y))  # unpenalized
    return {"slope": float(f.coef_[0, 0]), "intercept": float(f.intercept_[0]), "mean_score": float(np.mean(s)),
            "positive_rate": float(np.mean(y))}


def recall_threshold(y, s, recall: float) -> float:
    """Highest score at which recall reaches ``recall``."""
    y, s = np.asarray(y), np.asarray(s, float); o = np.argsort(-s)
    rec = np.cumsum(y[o]) / y.sum()
    return float(s[o][np.argmax(rec >= recall)])


def permutation_importance(predict: Callable[[pd.DataFrame], np.ndarray], table: pd.DataFrame, y,
                           groups: dict[str, Iterable[str]] = GROUPS, repeats: int = 20, seed: int = 0,
                           aggregate: Callable[[np.ndarray], np.ndarray] | None = None) -> tuple[np.ndarray, dict]:
    """Grouped permutation importance. Rows of each group's columns are permuted jointly; returns the unpermuted
    scores and ``{group: [scores per repeat]}``. ``aggregate`` maps row scores to pair scores when a pair has
    several rows (e.g. transcripts)."""
    agg = aggregate or (lambda p: p)
    rng = np.random.default_rng(seed)
    if any("DxA" in c for c in groups.values()) and "DxA" not in table and {"log10_distance", "Accessibility_Enh"} <= set(table):
        table = table.assign(DxA=table.log10_distance * table.Accessibility_Enh)   # permuted as a column, as in the manuscript
    full = agg(predict(table)); perm = {}
    for g, cols in groups.items():
        cols = [c for c in cols if c in table.columns]
        if not cols:
            continue
        reps = []
        for _ in range(repeats):
            t = table.copy(); ii = rng.permutation(len(t)); t[cols] = table[cols].to_numpy()[ii]
            reps.append(agg(predict(t)))
        perm[g] = reps
    return full, perm


def importance_table(y, clusters, full, perm, b: int = 500, seed: int = SEED, subset=None) -> pd.DataFrame:
    """Importance = AUPRC(unpermuted) - mean AUPRC over repeats, with a gene-clustered bootstrap CI
    (permuted scores held fixed). ``subset``: optional index array (e.g. one distance bin)."""
    y = np.asarray(y); idx = np.arange(len(y)) if subset is None else np.asarray(subset)
    bs = [idx[i] for i in cluster_bootstrap(np.asarray(clusters)[idx], b, seed, y[idx])]
    imp = lambda ii, reps: auprc(y[ii], full[ii]) - np.mean([auprc(y[ii], r[ii]) for r in reps])  # noqa: E731
    rows = []
    for g, reps in perm.items():
        v = [imp(ii, reps) for ii in bs]
        rows.append({"group": g, "importance": imp(idx, reps), "lo": np.percentile(v, 2.5), "hi": np.percentile(v, 97.5)})
    return pd.DataFrame(rows)
