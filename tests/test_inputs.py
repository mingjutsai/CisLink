"""Input handling: chromosome naming, expression formats, variants, preflight messages."""
import shutil

import numpy as np
import pandas as pd
import pytest
from test_features import BEDGRAPH, GTF_LINES, manual_sum

from cislink import features as cf
from cislink import inputs as ci

needs_bedtools = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")


def write_bedgraph(path, bare=False):
    path.write_text("".join(f"{c[3:] if bare else c}\t{s}\t{e}\t{v}\n" for c, s, e, v in BEDGRAPH))
    return str(path)


@pytest.fixture
def cfg(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text("".join(f"{c}\tT\t{f}\t{s}\t{e}\t.\t{st}\t.\t{i}\n" for c, f, s, e, st, i in GTF_LINES))
    pd.DataFrame({"transcript_id": ["T1.1", "T3.1", "T4.1"], "TPM": [10.0, 5.0, 2.0]}).to_csv(tmp_path / "tx.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": ["G1.1", "G2.1"], "TPM": [15.0, 2.0]}).to_csv(tmp_path / "g.tsv", sep="\t", index=False)
    (tmp_path / "cand.bed").write_text("chr1\t9000\t11000\t5\n")
    return {"gtf": str(gtf), "transcript_quantification": str(tmp_path / "tx.tsv"), "gene_quantification": str(tmp_path / "g.tsv"),
            "tracks": {"DNase": write_bedgraph(tmp_path / "d.bedGraph"), "H3K27ac": write_bedgraph(tmp_path / "k.bedGraph")},
            "candidates": str(tmp_path / "cand.bed"), "dir": tmp_path}


def test_naming_helpers():
    assert ci.strip_chr("chrX") == "X" and ci.format_chr("1", True) == "chr1" and ci.format_chr("chr1", False) == "1"
    assert ci.normalize_ensembl_id("ENSG00000123456.11") == "ENSG00000123456"


def test_expression_formats(tmp_path):
    rsem = tmp_path / "rsem.tsv"
    pd.DataFrame({"transcript_id": ["ENST1.2", "ENST2.1", "ENST2.1_PAR_Y"], "gene_id": ["G", "G", "G"], "TPM": [1.0, 2.0, 9.0]}).to_csv(rsem, sep="\t", index=False)
    two = tmp_path / "two.tsv"; two.write_text("ENST1.5\t3.0\nENST2\t4.0\n")
    hdr = tmp_path / "hdr.tsv"; hdr.write_text("id\tsample\nENST1\t5.0\n")
    assert ci.read_expression(rsem, "transcript").set_index("id").TPM.to_dict() == {"ENST1": 1.0, "ENST2": 2.0}
    assert ci.read_expression(two, "transcript").set_index("id").TPM.to_dict() == {"ENST1": 3.0, "ENST2": 4.0}
    assert ci.read_expression(hdr, "transcript").set_index("id").TPM.to_dict() == {"ENST1": 5.0}


def test_read_variants(tmp_path):
    a = tmp_path / "a.tsv"; a.write_text("SNPID\tchr\tpos\trs1\tchr1\t100\n".replace("\trs1", "\nrs1"))
    b = tmp_path / "b.tsv"; b.write_text("1\t100\trs1\n2\t200\n".replace("\n2\t200\n", "\n"))
    c = tmp_path / "c.tsv"; c.write_text("1\t100\n")
    assert ci.read_variants(a).iloc[0].tolist() == ["rs1", "chr1", 100]
    assert ci.read_variants(b).iloc[0].tolist() == ["rs1", "1", 100]
    assert ci.read_variants(c).iloc[0].tolist() == ["1:100", "1", 100]


@needs_bedtools
def test_bare_chromosome_tracks_give_the_same_signal(tmp_path):
    regions = pd.DataFrame({"chrom": ["chr1"] * 3, "start": [9750, 39750, 70000], "end": [10250, 40250, 70500]})
    s = cf.signal_sums(regions, {"DNase": write_bedgraph(tmp_path / "chr.bg"), "H3K27ac": write_bedgraph(tmp_path / "bare.bg", bare=True)})
    expected = [manual_sum(a, b) for a, b in zip(regions.start, regions.end)]
    assert np.allclose(s.DNase, expected) and np.allclose(s.H3K27ac, expected)


def test_variant_chromosomes_follow_the_annotation(cfg):
    tss = cf.protein_coding_tss(cf.read_gtf(cfg["gtf"])[0])
    r = cf.variant_pairs(pd.DataFrame({"SNPID": ["v"], "chr": ["1"], "pos": [40250]}), tss)
    assert len(r) == 3 and set(r.chrom) == {"chr1"}


@needs_bedtools
def test_preflight_accepts_and_rejects(cfg, tmp_path):
    info = ci.preflight(cfg, pd.DataFrame({"SNPID": ["v"], "chr": ["1"], "pos": [40250]}), log=lambda s: None)
    assert info["gtf_chr_prefix"] is True and info["track_chr_prefix"] == {"DNase": True, "H3K27ac": True}
    bad = dict(cfg, tracks={**cfg["tracks"], "H3K27ac": write_bedgraph(tmp_path / "bare.bg", bare=True)})
    with pytest.raises(ci.InputError, match="disagree on chromosome naming"):
        ci.preflight(bad, log=lambda s: None)
    with pytest.raises(ci.InputError, match="unknown track names"):
        ci.preflight(dict(cfg, tracks={**cfg["tracks"], "H3K27Ac": cfg["tracks"]["DNase"]}), log=lambda s: None)
    with pytest.raises(ci.InputError, match="accessibility track"):
        ci.preflight(dict(cfg, tracks={"H3K27ac": cfg["tracks"]["H3K27ac"]}), log=lambda s: None)
    (tmp_path / "none.tsv").write_text("ENST999.1\t1.0\n")
    with pytest.raises(ci.InputError, match="share no stable Ensembl IDs"):
        ci.preflight(dict(cfg, transcript_quantification=str(tmp_path / "none.tsv")), log=lambda s: None)
    with pytest.raises(ci.InputError, match="duplicated variant IDs"):
        ci.preflight(cfg, pd.DataFrame({"SNPID": ["v", "v"], "chr": ["1", "1"], "pos": [1, 2]}), log=lambda s: None)
