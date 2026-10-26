"""F1 features (cislink.roadmap) on a tiny synthetic Roadmap activity reference."""
import hashlib

import numpy as np
import pandas as pd
import pytest

from cislink import roadmap as rm

EIDS = ["E001", "E002", "E003", "E004", "E005"]
GTF = [
    ("chr1", "gene", 'gene_id "ENSG1.3"; gene_type "protein_coding"; gene_name "A";'),
    ("chr1", "transcript", 'gene_id "ENSG1.3"; transcript_id "ENST1.2"; gene_type "protein_coding"; gene_name "A";'),
    ("chr1", "gene", 'gene_id "ENSG2.1"; gene_type "lncRNA"; gene_name "B";'),
    ("chr1", "gene", 'gene_id "ENSG3.1"; gene_type "protein_coding"; gene_name "B";'),
]


@pytest.fixture
def ref(tmp_path):
    rng = np.random.default_rng(1)
    (tmp_path / "tracks").mkdir()
    man = []
    for mark in ("DNase", "H3K27ac"):
        for e in EIDS:
            np.savez_compressed(tmp_path / "tracks" / f"{e}-{mark}.npz", chr1=rng.uniform(0, 5, 400).astype(np.float32))
            man.append({"track": f"{e}-{mark}", "epigenome": e, "mark": mark, "bin_bp": 25})
    pd.DataFrame(man).to_csv(tmp_path / "manifest.tsv", sep="\t", index=False)
    pd.DataFrame(rng.uniform(0, 50, (2, 5)), index=["ENSG1", "ENSG3"], columns=EIDS).rename_axis("gene_id").to_csv(
        tmp_path / "rna_rpkm_55epigenomes.tsv.gz", sep="\t")
    files = ["manifest.tsv", "rna_rpkm_55epigenomes.tsv.gz"] + [f"tracks/{m['track']}.npz" for m in man]
    pd.DataFrame({"file": files, "bytes": [(tmp_path / f).stat().st_size for f in files],
                  "md5": [hashlib.md5((tmp_path / f).read_bytes()).hexdigest() for f in files]}).to_csv(tmp_path / "files_md5.tsv", sep="\t", index=False)
    return tmp_path


def test_reference_path_errors(monkeypatch, tmp_path):
    monkeypatch.delenv(rm.ENV, raising=False)
    with pytest.raises(rm.ReferenceError, match="reference download"):
        rm.reference_path(None)
    with pytest.raises(rm.ReferenceError, match="manifest"):
        rm.reference_path(str(tmp_path))


def test_signal_integral_with_partial_bins(ref):
    ra = rm.RoadmapActivity(str(ref))
    b = np.load(ref / "tracks" / "E001-DNase.npz")["chr1"].astype(float)
    s = ra.signals(["chr1:30:110"], "DNase").loc["chr1:30:110", "E001"]
    want = b[1] * 20 / 25 + b[2] + b[3] + b[4] * 10 / 25          # bins [25,50) [50,75) [75,100) [100,125)
    assert np.isclose(s, want)
    assert np.isclose(ra.signals(["chr1:30:40"], "DNase").iloc[0, 0], b[1] * 10 / 25)


def test_f1_is_spearman_and_missing_is_zero(ref):
    ra = rm.RoadmapActivity(str(ref))
    els = ["chr1:100:600", "chr1:1000:1500", np.nan, "chr1:100:600"]
    f = ra.f1_table(els, ["ENSG1", "ENSG3", "ENSG1", "ENSG9"])
    for mark, col in rm.MARKS:
        sig = ra.signals(["chr1:100:600"], mark).iloc[0]
        expr = np.log1p(ra.rpkm.loc["ENSG1", ra.eids[mark]])
        assert np.isclose(f.loc[0, col], pd.Series(sig.values).corr(pd.Series(expr.values), method="spearman"))
        assert f.loc[2, col] == 0 and f.loc[3, col] == 0          # unmapped element, unknown gene


def test_gene_ids(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text("".join(f"{c}\tT\t{f}\t1\t10\t.\t+\t.\t{i}\n" for c, f, i in GTF))
    ids = rm.gene_ids(gtf, ["ENST1.5", "ENST9.1", "ENST9.1"], ["A", "B", "Z"])
    assert ids[0] == "ENSG1" and ids[1] == "ENSG3" and pd.isna(ids[2])   # transcript first; symbol prefers protein-coding


def test_add_f1_hg19_and_verify(ref, tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text("".join(f"{c}\tT\t{f}\t1\t10\t.\t+\t.\t{i}\n" for c, f, i in GTF))
    t = pd.DataFrame({"Enh": ["chr1:100:600"], "TX": ["ENST1.2"], "gene": ["A"]})
    out = rm.add_f1(t, gtf, str(ref), genome="hg19", log=lambda s: None)
    want = rm.RoadmapActivity(str(ref)).f1_table(["chr1:100:600"], ["ENSG1"])
    assert np.allclose(out[["F1_dnase", "F1_h3k27ac"]].values, want.values)
    assert rm.verify(str(ref)) == 2 + 10
    (ref / "manifest.tsv").write_text("corrupt")
    with pytest.raises(rm.ReferenceError, match="corrupt"):
        rm.verify(str(ref))


def test_release_info():
    info = rm.release_info()
    assert info["directory"] == "roadmap_activity_v1" and len(info["archives"]) == 3
    assert all(len(a["md5"]) == 32 and a["bytes"] > 0 for a in info["archives"])
