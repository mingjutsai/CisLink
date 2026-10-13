import numpy as np
import pandas as pd

from cislink import evaluate as ev
from cislink import harmonize as hz


def test_reference_quantiles():
    q = hz.k562_chromatin_quantiles()
    assert len(q) == 1001 and set(hz.MARKS) <= set(q.columns)
    h = hz.k562_hic_quantiles()
    assert len(h) == 1001 and [f"s{i}" for i in range(16)] == list(h.columns[1:])


def test_k562_maps_onto_itself():
    q = hz.k562_chromatin_quantiles()
    t = pd.DataFrame({"H3K27ac_Enh": [0.0, 50.0, 500.0], "H3K27ac_Prom": [10.0, 100.0, 1000.0]})
    out = hz.harmonize_chromatin(t, q, enh_len=[500, 500, 500], prom_len=[500, 500, 500])
    assert np.allclose(out.values, t.values, rtol=1e-6, atol=1e-9)
    hq = hz.k562_hic_quantiles(); v = np.array([0.5, 1.0, np.nan, 3.0]); d = np.array([0, 20000, 50000, 200000])
    m = hz.harmonize_hic(v, d, hq)
    assert np.isnan(m[2]) and np.allclose(m[[0, 1, 3]], v[[0, 1, 3]], rtol=1e-6)


def test_quantile_mapping_is_monotone():
    rng = np.random.default_rng(0)
    cell = hz.chromatin_quantiles({"H3K27ac": rng.lognormal(0, 2, 5000)})
    f = hz.quantile_mapper(cell, hz.k562_chromatin_quantiles(), "H3K27ac")
    v = np.sort(rng.lognormal(0, 2, 100))
    assert np.all(np.diff(f(v)) >= 0)


def test_metrics():
    y = np.array([1, 0, 1, 0, 0, 1]); s = np.array([0.9, 0.1, 0.8, 0.3, 0.2, 0.7]); cl = list("abcdef")
    assert ev.auprc(y, s) == 1.0
    a, lo, hi = ev.auprc_ci(y, s, cl, b=50)
    assert lo <= a <= hi
    assert ev.recall_threshold(y, s, 0.5) == 0.8
    assert list(ev.distance_bins([5e3, 2e4, 1e5, 5e5])) == ["<10 kb", "10-50 kb", "50-200 kb", ">200 kb"]
