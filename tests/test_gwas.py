import numpy as np
import pandas as pd

from cislink import gwas


def silver():
    return pd.read_csv("data/gwas/silver_standard.tsv.gz", sep="\t").rename(
        columns={"Credible set": "locus", "Fine-mapping region": "region", "Silver gene": "target"})


def test_expected_rr():
    assert gwas.expected_rr(1, 1) == 1.0
    assert np.isclose(gwas.expected_rr(4, 1), (1 + 1 / 2 + 1 / 3 + 1 / 4) / 4)


def test_ties_split_and_no_call():
    t = pd.DataFrame({"locus": ["a"] * 3 + ["b"] * 2, "region": ["r1"] * 3 + ["r2"] * 2, "target": [True, False, False, False, True],
                      "s": [0.5, 0.5, 0.1, 0.0, 0.0], "near": [-10, -20, -30, -5, -50]})
    p = gwas.per_locus(t, "s")
    assert list(p.top1) == [0.5, 0.0] and list(p.called) == [True, False]
    assert gwas.per_locus(t, "near", always=True).top1.tolist() == [1.0, 0.0]


def test_manuscript_numbers():
    t = silver()
    s = gwas.summarize(t, {"CisLink-V2G": False, "Nearest gene": True, "L2G": False}).set_index("method")
    assert t.locus.nunique() == 568
    assert round(s.loc["CisLink-V2G", "top1"], 3) == 0.526 and round(s.loc["Nearest gene", "top1"], 3) == 0.518
    assert round(s.loc["L2G", "top1"], 3) == 0.451 and round(s.loc["CisLink-V2G", "MRR_cc"], 3) == 0.587
    lo = t[~t.locus.isin(set(t[t.target & t["L2G training gene"]].locus))]
    assert lo.locus.nunique() == 432
    assert round(gwas.locus_pr(gwas.locus_calls(lo, "CisLink-V2G"))["AP"], 3) == 0.287
