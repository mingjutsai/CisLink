"""Accessibility -log10 P tracks from BAM (cislink pvalue). Need samtools and bedtools; the MACS2 test also MACS2."""
import os
import shutil
import subprocess

import pandas as pd
import pytest

from cislink import pvalue as pv

needs_tools = pytest.mark.skipif(shutil.which("samtools") is None or shutil.which("bedtools") is None, reason="samtools/bedtools not installed")
MACS2 = os.environ.get("CISLINK_MACS2", ".pixi/envs/macs2/bin/macs2")
needs_macs2 = pytest.mark.skipif(shutil.which(MACS2) is None, reason="MACS2 not installed (pixi install -e macs2)")
HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr1_alt\tLN:5000\n"


def sam(tmp_path, reads, header=HEADER, name="r.bam"):
    """reads: (name, flag, pos 1-based, read length, mapq)."""
    body = "".join(f"{q}\t{f}\tchr1\t{p}\t{m}\t{n}M\t*\t0\t0\t{'A' * n}\t{'I' * n}\n" for q, f, p, n, m in reads)
    s = tmp_path / "r.sam"; s.write_text(header + body)
    b = tmp_path / name
    subprocess.run(f"samtools sort -o {b} {s} && samtools index {b}", shell=True, check=True)
    return b


def bed(path):
    return pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e", "n", "q", "strand"])


@needs_tools
def test_atac_tn5_shift(tmp_path):
    # plus-strand read 1001..1050 (0-based start 1000) -> start 1004; minus-strand read 2001..2050 (end 2050) -> end 2045
    b = sam(tmp_path, [("a", 0, 1001, 50, 60), ("b", 16, 2001, 50, 60)])
    n = pv.cut_sites(b, tmp_path / "c.bed", assay="atac")
    t = bed(tmp_path / "c.bed").sort_values("s")
    assert n == 2 and list(t.s) == [1004, 2000] and list(t.e) == [1050, 2045]


@needs_tools
def test_dnase_no_shift_and_filters(tmp_path):
    reads = [("a", 0, 1001, 36, 60), ("dup", 1024, 1001, 36, 60), ("low", 0, 5001, 36, 5), ("sec", 256, 7001, 36, 60)]
    b = sam(tmp_path, reads)
    assert pv.cut_sites(b, tmp_path / "a.bed", assay="dnase") == 2          # duplicate and secondary dropped
    assert list(bed(tmp_path / "a.bed").sort_values("s").s) == [1000, 5000]
    assert pv.cut_sites(b, tmp_path / "b.bed", assay="dnase", mapq=30) == 1


@needs_tools
def test_marks_duplicates_when_header_has_no_marking(tmp_path):
    b = sam(tmp_path, [("a", 0, 1001, 36, 60), ("b", 0, 1001, 36, 60), ("c", 0, 3001, 36, 60)])
    assert not pv.has_duplicate_marks(b)
    m = pv.mark_duplicates(b, tmp_path / "m.bam")
    assert pv.has_duplicate_marks(m) and pv.cut_sites(m, tmp_path / "c.bed", assay="dnase") == 2


@needs_tools
@needs_macs2
def test_pvalue_track_peaks_over_background(tmp_path):
    reads = [(f"p{i}", 0 if i % 2 else 16, 50001 + (i % 40), 50, 60) for i in range(200)]       # a pile at ~50 kb
    reads += [(f"b{i}", 0, 1001 + i * 997 % 98000, 50, 60) for i in range(300)]               # sparse background
    b = sam(tmp_path, reads)
    out = tmp_path / "acc.bedGraph"
    r = pv.pvalue_track(b, out, assay="atac", mark_dups="no", macs2=MACS2, log=lambda s: None)
    t = pd.read_csv(out, sep="\t", header=None, names=["c", "s", "e", "v"])
    assert r["cut_sites"] == 500 and set(t.c) == {"chr1"} and t.e.max() <= 100000
    peak = t[(t.s <= 50050) & (t.e > 50050)].v.max()
    assert peak > 10 and peak > 10 * t[(t.e < 40000)].v.median()


@needs_tools
def test_roadmap_dnase_tags(tmp_path):
    """Roadmap DNase tags: 5' 36 bp, MAPQ >= 30, duplicate-flagged and secondary reads dropped."""
    reads = [("p", 0, 1001, 50, 60), ("m", 16, 2001, 50, 60), ("low", 0, 5001, 50, 20), ("dup", 1024, 7001, 50, 60), ("sec", 256, 9001, 50, 60)]
    b = sam(tmp_path, reads)
    n = pv.roadmap_tags(b, tmp_path / "t.bed", pv.chrom_sizes(b))
    t = bed(tmp_path / "t.bed").sort_values("s")
    assert n == 2 and list(t.s) == [1000, 2014] and list(t.e) == [1036, 2050]   # plus: 5' start; minus: 5' end


@needs_tools
def test_roadmap_dnase_read1_only_and_subsample(tmp_path):
    reads = [("a", 65, 1001, 50, 60), ("a", 129, 1201, 50, 60), ("b", 65, 3001, 50, 60), ("b", 129, 3201, 50, 60), ("c", 0, 5001, 50, 60)]
    b = sam(tmp_path, reads)
    assert pv.roadmap_tags(b, tmp_path / "t.bed", pv.chrom_sizes(b)) == 2            # paired BAM: read 1 (flag 64) only
    assert list(bed(tmp_path / "t.bed").sort_values("s").s) == [1000, 3000]
    assert pv.roadmap_tags(b, tmp_path / "u.bed", pv.chrom_sizes(b), max_tags=1) == 1   # subsample


def test_dnase_needs_control(tmp_path):
    with pytest.raises(ValueError, match="control"):
        pv.pvalue_track(tmp_path / "x.bam", tmp_path / "x.bedGraph", assay="dnase")


@needs_tools
@needs_macs2
def test_replicates_are_pooled(tmp_path):
    reads = [(f"p{i}", 0 if i % 2 else 16, 50001 + (i % 40), 50, 60) for i in range(100)]
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    a, b = sam(tmp_path / "a", reads), sam(tmp_path / "b", reads)
    r = pv.pvalue_track([a, b], tmp_path / "acc.bedGraph", assay="atac", mark_dups="no", macs2=MACS2, log=lambda s: None)
    assert r["cut_sites"] == 200


def test_pool_and_cap(tmp_path):
    parts = []
    for i in range(2):
        p = tmp_path / f"t{i}.bed"; p.write_text("".join(f"chr1\t{j}\t{j + 36}\tN\t1000\t+\n" for j in range(i * 10, i * 10 + 10))); parts.append(str(p))
    assert pv._pool_and_cap(parts, tmp_path / "all.bed", None) == 20
    parts = []
    for i in range(2):
        p = tmp_path / f"u{i}.bed"; p.write_text("".join(f"chr1\t{j}\t{j + 36}\tN\t1000\t+\n" for j in range(i * 10, i * 10 + 10))); parts.append(str(p))
    assert pv._pool_and_cap(parts, tmp_path / "cap.bed", 15) == 15 and len((tmp_path / "cap.bed").read_text().splitlines()) == 15
