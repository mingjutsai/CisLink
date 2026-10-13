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
    assert len(load_model("e2g").features) == 37 and len(load_model("CisLink-V2G").features) == 13


def test_matches_sklearn_pipeline():
    """The JSON model equals the scikit-learn pipeline it was exported from."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import FunctionTransformer, StandardScaler

    from cislink.model import slog
    for m in MODELS:
        d = json.load(resources.files("cislink.models").joinpath(f"{m}.json").open())
        pipe = make_pipeline(FunctionTransformer(slog), StandardScaler(), LogisticRegression())
        x = pd.read_csv("examples/example_features.tsv", sep="\t")[d["features"]].fillna(0.0)
        pipe.fit(np.vstack([x.values, x.values]), np.r_[np.zeros(len(x)), np.ones(len(x))])  # set shapes, then overwrite
        sc, lr = pipe.named_steps["standardscaler"], pipe.named_steps["logisticregression"]
        sc.mean_, sc.scale_ = np.array(d["scaler_mean"]), np.array(d["scaler_scale"])
        lr.coef_, lr.intercept_ = np.array([d["coef"]]), np.array([d["intercept"]])
        assert np.allclose(pipe.predict_proba(x.values)[:, 1], load_model(m).predict(x), atol=1e-12)


@pytest.mark.parametrize("cell", ["K562", "HCT116", "GM12878"])
def test_reproduces_reference_scores(cell):
    p = pd.read_csv(f"{DATA}/{cell}_heldout.pairs.tsv.gz", sep="\t")
    f = pd.read_csv(f"{DATA}/{cell}_heldout.features.tsv.gz", sep="\t")
    for m, ref in (("e2g", "reference CisLink-E2G"), ("v2g", "reference CisLink-V2G")):
        s = pd.Series(load_model(m).predict(f), index=f.ID).groupby(level=0).max()
        assert np.allclose(p.ID.map(s).values, p[ref].values, atol=1e-5)


def test_missing_feature_raises():
    x = pd.read_csv("examples/example_features.tsv", sep="\t").drop(columns=["H3K27ac_Enh"])
    with pytest.raises(KeyError):
        load_model("e2g").predict(x)
    load_model("v2g").predict(x)  # V2G does not use histone marks


def test_missing_hic_is_zero():
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    y = x.copy(); y["HiC_Contact"] = np.nan; z = x.copy(); z["HiC_Contact"] = 0.0
    assert np.allclose(load_model("e2g").predict(y), load_model("e2g").predict(z))


def test_contributions_sum_to_logit():
    x = pd.read_csv("examples/example_features.tsv", sep="\t"); m = load_model("e2g")
    assert np.allclose(m.contributions(x).sum(axis=1) + m.intercept, m.decision_function(x))


def test_score_and_cli(tmp_path):
    x = pd.read_csv("examples/example_features.tsv", sep="\t")
    s = score(x)
    assert {"CisLink-E2G", "CisLink-V2G"} <= set(s.columns) and s["CisLink-E2G"].between(0, 1).all()
    out = tmp_path / "s.tsv"
    assert main(["score", "examples/example_features.tsv", "-o", str(out), "--keep", "ID", "--threshold", "70"]) == 0
    t = pd.read_csv(out, sep="\t")
    assert list(t.columns[:3]) == ["ID", "CisLink-E2G", "CisLink-E2G >= 70% recall threshold"]
    assert np.allclose(t["CisLink-E2G"], s["CisLink-E2G"], atol=1e-6)
