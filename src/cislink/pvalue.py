"""Accessibility -log10 P track from a BAM file (``cislink pvalue``).

The models read accessibility as a MACS2 -log10 Poisson P track, the format of the training data. A BAM file is turned
into that track as the ENCODE ATAC-seq pipeline does it:

1. Reads: mapped, primary, not QC-failed, not duplicates (``-F 3844``), optionally ``MAPQ >= --mapq``. When the BAM has
   no duplicate marks (no ``samtools markdup`` / Picard ``MarkDuplicates`` in its header), duplicates are marked
   first (``samtools collate | fixmate -m | sort | markdup``; paired- or single-end).
2. Cut sites: ATAC-seq, the Tn5 insertion (5' end shifted +4 on the plus strand, -5 on the minus strand); DNase-seq,
   the 5' end. Each mate of a pair is a cut site.
3. ``macs2 callpeak -f BED --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01`` without a control,
   then ``macs2 bdgcmp -m ppois -S <cut sites / 1e6>``. The result is clipped to the BAM's chromosomes and sorted
   (bedGraph; also a bigWig when the output ends in .bw/.bigWig).

On osteoblast ATAC-seq the track agrees with the ENCODE ATAC pipeline's own p-value track (Spearman 0.996 over the
candidate elements). Needs ``samtools`` and ``bedtools`` (``pipeline`` pixi environment) and MACS2 2.2.7.1 (``macs2``
pixi environment; pass its path with ``--macs2``).
"""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile

ASSAYS = {"atac": (4, -5), "dnase": (0, 0)}   # 5'-end shift on the plus / minus strand
READ_FLAGS_EXCLUDED = 3844   # unmapped, secondary, QC fail, duplicate, supplementary
MACS2_CALLPEAK = "--nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01"


def _run(cmd: str, **kw) -> str:
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True, **kw).stdout


def has_duplicate_marks(bam, samtools: str = "samtools") -> bool:
    """True when the BAM header records a duplicate-marking step."""
    head = _run(f"{shlex.quote(samtools)} view -H {shlex.quote(str(bam))}")
    return any(line.startswith("@PG") and any(k in line for k in ("markdup", "MarkDuplicates", "rmdup")) for line in head.splitlines())


def is_paired(bam, samtools: str = "samtools", n: int = 10000) -> bool:
    flags = _run(f"{shlex.quote(samtools)} view {shlex.quote(str(bam))} | head -n {n} | cut -f2")
    return any(int(f) & 1 for f in flags.split())


def chrom_sizes(bam, samtools: str = "samtools") -> dict[str, int]:
    """Chromosomes of the BAM header without alternate, unplaced or decoy contigs (names without '_' or 'EBV')."""
    out = {}
    for line in _run(f"{shlex.quote(samtools)} view -H {shlex.quote(str(bam))}").splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            if "_" not in f["SN"] and "EBV" not in f["SN"] and f["SN"] not in ("chrM", "MT", "M"):
                out[f["SN"]] = int(f["LN"])
    return out


def mark_duplicates(bam, out_bam, samtools: str = "samtools", threads: int = 4, tmpdir=None) -> str:
    s, t = shlex.quote(samtools), int(threads)
    tmp = tempfile.mkdtemp(dir=tmpdir)
    try:
        pre = (f"{s} collate -@ {t} -O -u {shlex.quote(str(bam))} {tmp}/collate | {s} fixmate -@ {t} -m -u - - | "
               if is_paired(bam, samtools) else f"{s} view -@ {t} -u {shlex.quote(str(bam))} | ")
        _run(pre + f"{s} sort -@ {t} -u -T {tmp}/sort - | {s} markdup -@ {t} - {shlex.quote(str(out_bam))}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return str(out_bam)


def cut_sites(bam, out, assay: str = "atac", mapq: int = 0, samtools: str = "samtools", bedtools: str = "bedtools",
              threads: int = 4) -> int:
    """Write the cut sites of the filtered reads as a 6-column BED of the shifted reads (MACS2 input); returns their number."""
    plus, minus = ASSAYS[assay]
    view = f"{shlex.quote(samtools)} view -@ {int(threads)} -b -F {READ_FLAGS_EXCLUDED} -q {int(mapq)} {shlex.quote(str(bam))}"
    _run(f"{view} | {shlex.quote(bedtools)} bamtobed -i stdin"
         f" | awk -v OFS='\\t' '{{if ($6 == \"+\") $2 = $2 + {plus}; else $3 = $3 + {minus}; if ($3 > $2) print $1, $2, $3, \"N\", 1000, $6}}'"
         f" > {shlex.quote(str(out))}")
    return int(_run(f"wc -l < {shlex.quote(str(out))}").strip())


def pvalue_track(bam, output, assay: str = "atac", mapq: int = 0, mark_dups: str = "auto", macs2: str = "macs2",
                 samtools: str = "samtools", bedtools: str = "bedtools", threads: int = 4, tmpdir=None, log=print) -> dict:
    """Write the accessibility -log10 P track (bedGraph; also a bigWig when ``output`` ends in .bw/.bigWig) of a BAM
    file. ``mark_dups``: 'auto' (mark when the header shows no marking step), 'yes' or 'no'. Returns a summary."""
    assay = assay.lower()
    if assay not in ASSAYS:
        raise ValueError(f"assay must be one of {list(ASSAYS)}")
    for tool in (samtools, bedtools, macs2):
        if not shutil.which(tool):
            raise RuntimeError(f"{tool} not found (MACS2: pixi install -e macs2, then --macs2 .pixi/envs/macs2/bin/macs2)")
    out = str(output); bigwig = out.lower().endswith((".bw", ".bigwig"))
    bedgraph = (out.rsplit(".", 1)[0] + ".bedGraph") if bigwig else out
    tmp = tempfile.mkdtemp(dir=tmpdir or os.path.dirname(os.path.abspath(out)))
    try:
        src = str(bam)
        if mark_dups == "yes" or (mark_dups == "auto" and not has_duplicate_marks(bam, samtools)):
            log("marking duplicates (samtools fixmate/markdup)")
            src = mark_duplicates(bam, os.path.join(tmp, "markdup.bam"), samtools, threads, tmp)
        sizes = chrom_sizes(src, samtools)
        tags = os.path.join(tmp, "cuts.bed")
        n = cut_sites(src, tags, assay, mapq, samtools, bedtools, threads)
        if n == 0:
            raise RuntimeError("no reads pass the filters")
        log(f"{n} cut sites ({assay}); MACS2")
        m, t = shlex.quote(macs2), shlex.quote(tmp)
        _run(f"{m} callpeak -t {shlex.quote(tags)} -f BED -n acc -g hs --outdir {t} {MACS2_CALLPEAK}")
        _run(f"{m} bdgcmp -t {t}/acc_treat_pileup.bdg -c {t}/acc_control_lambda.bdg -m ppois -S {n / 1e6!r} --o-prefix {t}/acc")
        cs = os.path.join(tmp, "chrom.sizes")
        with open(cs, "w") as fh:
            fh.writelines(f"{c}\t{size}\n" for c, size in sizes.items())
        clip = "'NR == FNR {L[$1] = $2; next} ($1 in L) {if ($3 > L[$1]) $3 = L[$1]; if ($3 > $2) print}'"
        _run(f"awk -v OFS='\\t' {clip} {shlex.quote(cs)} {t}/acc_ppois.bdg"
             f" | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T {t} > {shlex.quote(bedgraph)}")
        log(f"-log10 P track -> {bedgraph}")
        if bigwig:
            import pyBigWig
            bw = pyBigWig.open(out, "w"); bw.addHeader(list(sizes.items()))
            for chunk in _bedgraph_chunks(bedgraph):
                bw.addEntries(*chunk)
            bw.close(); log(f"bigWig -> {out}")
        return {"cut_sites": n, "assay": assay, "bedgraph": bedgraph, "bigwig": out if bigwig else None}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _bedgraph_chunks(path, size: int = 1_000_000):
    c, s, e, v = [], [], [], []
    with open(path) as fh:
        for line in fh:
            a = line.split("\t")
            c.append(a[0]); s.append(int(a[1])); e.append(int(a[2])); v.append(float(a[3]))
            if len(c) >= size:
                yield c, s, e, v
                c, s, e, v = [], [], [], []
    if c:
        yield c, s, e, v
