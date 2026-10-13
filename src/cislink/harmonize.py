"""Put a new cell type's inputs on the K562 scale the models were trained on.

* Chromatin tracks (accessibility, histone marks, CTCF): the per-bp mean signal over the cell's
  genome-wide candidate elements is summarized as 1,001 quantiles; each element or promoter sum is
  converted to a per-bp mean, mapped quantile-to-quantile onto K562, and multiplied by the region
  length again. A region that is closed in the cell stays low; only the assay scale is aligned.
* Expression: TPM recomputed over protein-coding genes and transcripts (multiply by 1e6 / the
  protein-coding TPM sum).
* Hi-C: observed/expected contact mapped onto the K562 5-kb distribution within 16 distance strata.

The K562 references ship with the package (``reference/``).
"""
from __future__ import annotations

from importlib import resources

import numpy as np
import pandas as pd

MARKS = ("Accessibility", "H4K20me1", "H3K9me3", "H3K9ac", "H3K79me2", "H3K4me3", "H3K4me2", "H3K4me1", "H3K36me3",
         "H3K27me3", "H3K27ac", "H2A.Z", "CTCF")
Q = np.linspace(0, 1, 1001)
HIC_EDGES = np.array([0, 5e3, 10e3, 15e3, 20e3, 30e3, 40e3, 50e3, 75e3, 100e3, 150e3, 200e3, 300e3, 400e3, 500e3, 750e3, 1e6 + 1])


def _ref(name: str) -> pd.DataFrame:
    with resources.files("cislink.reference").joinpath(name).open() as fh:
        return pd.read_csv(fh, sep="\t")


def compat(q: pd.DataFrame) -> pd.DataFrame:
    """Quantile tables written before accessibility was renamed ('DNase') still work."""
    return q.rename(columns={"DNase": "Accessibility"}) if "DNase" in q and "Accessibility" not in q else q


def k562_chromatin_quantiles() -> pd.DataFrame:
    """1,001 quantiles of the per-bp mean signal of each track over K562 candidate elements."""
    return _ref("K562_chromatin_quantiles.tsv")


def k562_hic_quantiles() -> pd.DataFrame:
    """1,001 quantiles of K562 5-kb Hi-C observed/expected per distance stratum (columns s0..s15)."""
    return _ref("K562_hic_oe_5kb_quantiles.tsv")


def chromatin_quantiles(per_bp_means: dict[str, np.ndarray]) -> pd.DataFrame:
    """Quantiles of a cell's tracks: ``{track: per-bp mean signal over the cell's candidate elements}``.
    Tracks the cell lacks are left out (and are then not mapped)."""
    out = {"q": Q}
    for m, v in per_bp_means.items():
        v = np.asarray(v, float)
        out[m] = np.quantile(v[~np.isnan(v)], Q)
    return pd.DataFrame(out)


def quantile_mapper(cell_q: pd.DataFrame, ref_q: pd.DataFrame, track: str):
    """Function mapping a cell's per-bp mean signal of ``track`` onto the reference scale."""
    x, y = cell_q[track].to_numpy(float), ref_q[track].to_numpy(float)
    x, idx = np.unique(x, return_index=True)  # ties at the low end (zeros)
    return lambda v: np.interp(v, x, y[idx])


def harmonize_chromatin(table: pd.DataFrame, cell_q: pd.DataFrame, enh_len, prom_len, ref_q: pd.DataFrame | None = None) -> pd.DataFrame:
    """Map every ``{track}_Enh`` / ``{track}_Prom`` sum in ``table`` onto K562. ``enh_len`` and ``prom_len``
    are the region lengths in bp (array-like, row-aligned)."""
    ref_q = k562_chromatin_quantiles() if ref_q is None else compat(ref_q)
    cell_q = compat(cell_q)
    out = table.copy()
    we = np.clip(np.asarray(enh_len, float), 1, None); wp = np.clip(np.asarray(prom_len, float), 1, None)
    for m in MARKS:
        if m not in cell_q or cell_q[m].isna().all():
            continue
        f = quantile_mapper(cell_q, ref_q, m)
        for part, w in (("Enh", we), ("Prom", wp)):
            col = f"{m}_{part}"
            if col in out:
                out[col] = f(out[col].to_numpy(float) / w) * w
    return out


def expression_factors(gene_tpm_protein_coding_sum: float, transcript_tpm_protein_coding_sum: float) -> dict:
    """Scale factors that make expression protein-coding TPM."""
    return {"gene_factor": 1e6 / gene_tpm_protein_coding_sum, "transcript_factor": 1e6 / transcript_tpm_protein_coding_sum}


def harmonize_expression(table: pd.DataFrame, factors: dict) -> pd.DataFrame:
    out = table.copy()
    out["Gene_expression"] = out["Gene_expression"] * factors["gene_factor"]
    out["Tx_expression"] = out["Tx_expression"] * factors["transcript_factor"]
    return out


def hic_stratum(bin_distance) -> np.ndarray:
    d = np.asarray(bin_distance, float)
    return np.clip(np.searchsorted(HIC_EDGES, d, side="right") - 1, 0, len(HIC_EDGES) - 2)


def hic_quantiles(oe: np.ndarray, bin_distance: np.ndarray) -> pd.DataFrame:
    """A cell's Hi-C O/E quantiles per distance stratum, from sampled bin pairs (values and bin distances)."""
    oe, st = np.asarray(oe, float), hic_stratum(bin_distance)
    out = {"q": Q}
    for s in range(len(HIC_EDGES) - 1):
        v = oe[(st == s) & ~np.isnan(oe)]
        out[f"s{s}"] = np.quantile(v, Q) if len(v) else np.full(len(Q), np.nan)
    return pd.DataFrame(out)


def harmonize_hic(oe, bin_distance, cell_q: pd.DataFrame, ref_q: pd.DataFrame | None = None) -> np.ndarray:
    """Map Hi-C O/E values (element bin to TSS bin) onto the K562 5-kb scale, stratum by stratum.
    Missing values stay missing (the models treat them as 0)."""
    ref_q = k562_hic_quantiles() if ref_q is None else ref_q
    v = np.asarray(oe, float); out = np.full(len(v), np.nan); st = hic_stratum(bin_distance)
    for s in np.unique(st):
        k = (st == s) & ~np.isnan(v)
        c, r = cell_q[f"s{s}"].to_numpy(float), ref_q[f"s{s}"].to_numpy(float)
        if k.any() and not np.isnan(c).all():
            out[k] = np.interp(v[k], c, r)
    return out
