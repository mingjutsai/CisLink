"""Feature pipeline on a tiny synthetic genome. bedGraph tests need bedtools, bigWig tests need pyBigWig."""
import shutil

import numpy as np
import pandas as pd
import pytest

from cislink import features as cf
from cislink import harmonize as hz
from cislink import load_model

GTF_LINES = [
    # gene A (+) with two transcripts sharing a TSS and one with another TSS; gene B (-); gene C (+) non-coding
    ("chr1", "gene", 10000, 20000, "+", 'gene_id "G1.1"; gene_type "protein_coding"; gene_name "A";'),
    ("chr1", "transcript", 10000, 20000, "+", 'gene_id "G1.1"; transcript_id "T1.1"; gene_type "protein_coding"; gene_name "A"; transcript_type "protein_coding";'),
    ("chr1", "transcript", 10000, 18000, "+", 'gene_id "G1.1"; transcript_id "T2.1"; gene_type "protein_coding"; gene_name "A"; transcript_type "retained_intron";'),
    ("chr1", "transcript", 12000, 20000, "+", 'gene_id "G1.1"; transcript_id "T3.1"; gene_type "protein_coding"; gene_name "A"; transcript_type "protein_coding";'),
    ("chr1", "gene", 30000, 60000, "-", 'gene_id "G2.1"; gene_type "protein_coding"; gene_name "B";'),
    ("chr1", "transcript", 30000, 60000, "-", 'gene_id "G2.1"; transcript_id "T4.1"; gene_type "protein_coding"; gene_name "B"; transcript_type "protein_coding";'),
    ("chr1", "gene", 80000, 90000, "+", 'gene_id "G3.1"; gene_type "lncRNA"; gene_name "C";'),
    ("chr1", "transcript", 80000, 90000, "+", 'gene_id "G3.1"; transcript_id "T5.1"; gene_type "lncRNA"; gene_name "C"; transcript_type "lncRNA";'),
]
BEDGRAPH = [("chr1", 0, 9800, 0.0), ("chr1", 9800, 10100, 2.0), ("chr1", 10100, 40000, 0.5), ("chr1", 40000, 40500, 4.0), ("chr1", 40500, 100000, 0.0)]


@pytest.fixture
def inputs(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text("".join(f"{c}\tT\t{f}\t{s}\t{e}\t.\t{st}\t.\t{i}\n" for c, f, s, e, st, i in GTF_LINES))
    bg = tmp_path / "s.bedGraph"
    bg.write_text("".join(f"{c}\t{s}\t{e}\t{v}\n" for c, s, e, v in BEDGRAPH))
    tq = tmp_path / "tx.tsv"
    pd.DataFrame({"transcript_id": ["T1.1", "T3.1", "T4.1", "T2.1"], "TPM": [10.0, 5.0, 2.0, 1.0]}).to_csv(tq, sep="\t", index=False)
    gq = tmp_path / "gene.tsv"
    pd.DataFrame({"gene_id": ["G1.1", "G2.1", "G3.1"], "TPM": [15.0, 2.0, 7.0]}).to_csv(gq, sep="\t", index=False)
    return {"gtf": gtf, "bedgraph": bg, "tq": tq, "gq": gq, "dir": tmp_path}


def manual_sum(s, e):
    return sum(v * max(0, min(e, b) - max(s, a)) for _, a, b, v in BEDGRAPH)


def test_gtf_and_tss(inputs):
    tx, genes = cf.read_gtf(inputs["gtf"])
    assert list(tx.TX) == ["T1.1", "T2.1", "T3.1", "T4.1", "T5.1"]
    tss = cf.protein_coding_tss(tx)
    assert list(tss.TX) == ["T1.1", "T3.1", "T4.1"]   # T2 shares T1's TSS; T5 is not protein-coding
    assert list(tss.tss) == [10000, 12000, 60000]     # minus-strand TSS is the end


def test_variant_pairs(inputs):
    tss = cf.protein_coding_tss(cf.read_gtf(inputs["gtf"])[0])
    r = cf.variant_pairs(pd.DataFrame({"SNPID": ["v1"], "chr": ["chr1"], "pos": [40250]}), tss, max_distance=30_000)
    assert list(r.TX) == ["T3.1", "T4.1"]
    assert r.ID.iloc[0] == "v1:A|T3.1" and (r.ee - r.es == 500).all() and (r.pe - r.ps == 500).all()


def test_distance_features(inputs):
    tx, genes = cf.read_gtf(inputs["gtf"])
    r = cf.variant_pairs(pd.DataFrame({"SNPID": ["v1"], "chr": ["chr1"], "pos": [40250]}), cf.protein_coding_tss(tx))
    d = cf.distance_features(r, genes).set_index(r.TX)
    assert d.loc["T4.1", "distance"] == 19750 and d.loc["T4.1", "gene_rank"] == 1 and d.loc["T4.1", "is_nearest_gene"] == 1
    assert d.loc["T1.1", "gene_rank"] == 2 and d.loc["T1.1", "n_genes_between"] == 0
    assert np.isclose(d.loc["T4.1", "powerlaw_contact"], (19750 / 5000) ** -1.024)


def test_expression(inputs):
    tx, genes = cf.read_gtf(inputs["gtf"])
    assert list(cf.transcript_expression(["T1.1", "T9.1"], inputs["tq"])) == [10.0, 0.0]
    assert list(cf.gene_expression(["A", "B", "Z"], genes, inputs["gq"])) == [15.0, 2.0, 0.0]
    f = cf.expression_factors(tx, genes, inputs["tq"], inputs["gq"])
    assert np.isclose(f["gene_factor"], 1e6 / 17.0) and np.isclose(f["transcript_factor"], 1e6 / 17.0)   # T2 is not protein-coding


REGIONS = pd.DataFrame({"chrom": ["chr1"] * 3, "start": [9750, 39750, 70000], "end": [10250, 40250, 70500]})


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")
def test_signal_sums_bedgraph(inputs):
    s = cf.signal_sums(REGIONS, {"accessibility": str(inputs["bedgraph"])})
    assert np.allclose(s.Accessibility, [manual_sum(a, b) for a, b in zip(REGIONS.start, REGIONS.end)])
    assert s.H3K27ac.isna().all()


def test_signal_sums_bigwig(inputs):
    pbw = pytest.importorskip("pyBigWig")
    path = str(inputs["dir"] / "s.bw")
    bw = pbw.open(path, "w"); bw.addHeader([("chr1", 100000)])
    bw.addEntries(["chr1"] * len(BEDGRAPH), [a for _, a, _, _ in BEDGRAPH], ends=[b for _, _, b, _ in BEDGRAPH], values=[v for *_, v in BEDGRAPH])
    bw.close()
    s = cf.signal_sums(REGIONS, {"DNase": path})   # earlier name still accepted
    assert np.allclose(s.Accessibility, [manual_sum(a, b) for a, b in zip(REGIONS.start, REGIONS.end)])


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")
def test_build_without_hic(inputs):
    tx, _ = cf.read_gtf(inputs["gtf"])
    rows = cf.variant_pairs(pd.DataFrame({"SNPID": ["v1", "v2"], "chr": ["chr1", "chr1"], "pos": [40250, 11000]}), cf.protein_coding_tss(tx))
    tracks = {m: str(inputs["bedgraph"]) for m in hz.MARKS}
    out = cf.build_features(rows, gtf=inputs["gtf"], tracks=tracks, transcript_quantification=inputs["tq"], gene_quantification=inputs["gq"],
                            candidates_bed=None, cell_quantiles=hz.k562_chromatin_quantiles(), hic=None, megamap=None, log=lambda s: None)
    for m in ("v2g", "e2g"):
        p = load_model(m).predict(out)
        assert len(p) == len(rows) and np.all((p >= 0) & (p <= 1))


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")
def test_signal_sums_parallel_tracks_agree(inputs):
    """Regression: tracks summed in parallel must not share temporary files (each must equal the serial result)."""
    regions = pd.DataFrame({"chrom": ["chr1"] * 200, "start": np.arange(200) * 400, "end": np.arange(200) * 400 + 500})
    expected = [manual_sum(a, b) for a, b in zip(regions.start, regions.end)]
    for _ in range(3):
        s = cf.signal_sums(regions, {m: str(inputs["bedgraph"]) for m in hz.MARKS}, threads=len(hz.MARKS))
        for m in hz.MARKS:
            assert np.allclose(s[m], expected), m
