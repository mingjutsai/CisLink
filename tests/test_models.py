import json
from importlib import resources

import numpy as np
import pandas as pd
import pytest

from cislink import MODELS, load_model, score
from cislink.cli import main

DATA = "data/crispr"


def test_models_load():
    for m in MODELS:
        model = load_model(m)
        assert len(model.features) == len(model.coef) == len(model.scaler_mean)
    assert len(load_model("extended").features) == 40 and len(load_model("CisLink-Base").features) == 16


def test_matches_sklearn_pipeline():
    """The JSON model equals the scikit-learn pipeline it was exported from."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import FunctionTransformer, StandardScaler

    from cislink.model import slog
    for m in MODELS:
        d = json.load(resources.files("cislink.models").joinpath(f"{m}.json").open())
        pipe = make_pipeline(FunctionTransformer(slog), StandardScaler(), LogisticRegression())
        x = pd.read_csv("examples/example_features.tsv", sep="\t")
        x = x.assign(DxA=x.log10_distance * x.Accessibility_Enh)[d["features"]].fillna(0.0)
        pipe.fit(np.vstack([x.values, x.values]), np.r_[np.zeros(len(x)), np.ones(len(x))])  # set shapes, then overwrite
        sc, lr = pipe.named_steps["standardscaler"], pipe.named_steps["logisticregression"]
        sc.mean_, sc.scale_ = np.array(d["scaler_mean"]), np.array(d["scaler_scale"])
        lr.coef_, lr.intercept_ = np.array([d["coef"]]), np.array([d["intercept"]])
        assert np.allclose(pipe.predict_proba(x.values)[:, 1], load_model(m).predict(x), atol=1e-12)


@pytest.mark.parametrize("cell", ["K562", "HCT116", "GM12878", "WTC11"])
def test_reproduces_reference_scores(cell):
    p = pd.read_csv(f"{DATA}/{cell}_heldout.pairs.tsv.gz", sep="\t")
    f = pd.read_csv(f"{DATA}/{cell}_heldout.features.tsv.gz", sep="\t")
    for m, ref in (("extended", "reference CisLink-Extended"), ("base", "reference CisLink-Base")):
        if p[ref].isna().all():   # WTC11: no histone marks, CisLink-Extended not scored
            continue
        s = pd.Series(load_model(m).predict(f), index=f.ID).groupby(level=0).max()
        assert np.allclose(p.ID.map(s).values, p[ref].values, atol=1e-5)


def test_crispr_manuscript_numbers():
    """Table 1a AUPRC of the cells outside K562, from the shipped feature rows."""
    from sklearn.metrics import average_precision_score as ap

    from cislink import datasets
    want = {"HCT116": (0.750, 0.513), "GM12878": (0.913, 0.880), "WTC11": (None, 0.429)}
    for cell, (ext, bas) in want.items():
        p = datasets.crispr_pairs(cell, path=DATA)
        assert round(ap(p.label, p["CisLink-Base"]), 3) == bas
        if ext is None:
            assert p["CisLink-Extended"].isna().all()
        else:
            assert round(ap(p.label, p["CisLink-Extended"]), 3) == ext


def test_missing_feature_raises():
    x = pd.read_csv("examples/example_features.tsv", sep="\t").drop(columns=["H3K27ac_Enh"])
    with pytest.raises(KeyError):
        load_model("extended").predict(x)
    load_model("base").predict(x)  # CisLink-Base does not use histone marks


def test_missing_hic_is_zero():
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    y = x.copy(); y["HiC_Contact"] = np.nan; z = x.copy(); z["HiC_Contact"] = 0.0
    assert np.allclose(load_model("extended").predict(y), load_model("extended").predict(z))


def test_contributions_sum_to_logit():
    x = pd.read_csv("examples/example_features.tsv", sep="\t"); m = load_model("extended")
    assert np.allclose(m.contributions(x).sum(axis=1) + m.intercept, m.decision_function(x))


def test_score_and_cli(tmp_path):
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    s = score(x)
    assert {"CisLink-Extended", "CisLink-Base"} <= set(s.columns) and s["CisLink-Extended"].between(0, 1).all()
    out = tmp_path / "s.tsv"
    assert main(["score", "examples/example_features.tsv", "-o", str(out), "--keep", "ID", "--threshold", "70"]) == 0
    t = pd.read_csv(out, sep="\t")
    assert list(t.columns[:3]) == ["ID", "CisLink-Extended", "CisLink-Extended >= 70% recall threshold"]
    assert np.allclose(t["CisLink-Extended"], s["CisLink-Extended"], atol=1e-6)


def test_released_version_and_thresholds():
    """CisLink v1.3 models (renamed in 1.4.0, unchanged): the frozen models (source hashes) and the 70% / 50% recall
    thresholds from cross-validated training predictions."""
    e, v = load_model("extended"), load_model("base")
    assert e.meta["version"] == v.meta["version"] == "1.3.0"
    assert e.meta["source_sha256"].startswith("2a1fad0e") and v.meta["source_sha256"].startswith("95977c92")
    assert (e.threshold(70), e.threshold(50), v.threshold(70), v.threshold(50)) == (0.269, 0.565, 0.23, 0.476)


def test_k562_heldout_auprc():
    from sklearn.metrics import average_precision_score as ap

    from cislink import datasets
    p = datasets.crispr_pairs("K562", path=DATA)
    assert (round(ap(p.label, p["CisLink-Extended"]), 3), round(ap(p.label, p["CisLink-Base"]), 3)) == (0.700, 0.603)


def test_dxa_is_derived_and_f1_missing_is_zero():
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    for m in MODELS:
        model = load_model(m)
        with_dxa = x.assign(DxA=x.log10_distance * x.Accessibility_Enh)
        assert np.array_equal(model.predict(x), model.predict(with_dxa))
        y = x.copy(); y["F1_dnase"] = np.nan; z = x.copy(); z["F1_dnase"] = 0.0
        assert np.allclose(model.predict(y), model.predict(z))
    with pytest.raises(KeyError, match="add-f1"):
        load_model("base").predict(x.drop(columns=["F1_dnase", "F1_h3k27ac"]))


def test_legacy_accessibility_columns():
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    old = x.rename(columns={"Accessibility_Enh": "DNase_Enh", "Accessibility_Prom": "DNase_Prom"})
    for m in MODELS:
        assert np.array_equal(load_model(m).predict(x), load_model(m).predict(old))

