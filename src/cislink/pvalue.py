"""Accessibility -log10 P track from a BAM file (``cislink pvalue``).

The models read accessibility as a MACS2 -log10 Poisson P track, the format of the training data (Roadmap consolidated
DNase-seq). Two procedures, one per assay:

DNase-seq (``--assay dnase``): the Roadmap Epigenomics procedure that produced the training track, reconstructed and
validated on K562 (element-level Spearman 0.998 with the Roadmap E123 track):

1. Tags: mapped, primary, not QC-failed, not duplicate-flagged reads with MAPQ >= 30; read 1 only for paired-end data;
   each read truncated to its 5' 36 bp (Roadmap's uniform read length); at most 50 million tags (random subsample,
   fixed seed).
   Several BAMs (biological replicates) are pooled before the cap, as in Roadmap's consolidated tracks.
2. Control: Roadmap's simulated uniform DNase control (E123-DNaseControl, lifted to GRCh38; ``cislink roadmap-control``
   downloads and lifts it).
3. ``macs2 callpeak -t tags -c control -f BED -g hs -p 0.01 --nomodel --extsize 146 -B --SPMR`` (MACS2 2.2.7.1; equal
   to Roadmap's MACS2 2.0.10 ``--shiftsize 73``; duplicates at MACS2's default, one per position and strand), then
   ``macs2 bdgcmp -m ppois -S <min(tags, control tags) / 1e6>``.

ATAC-seq (``--assay atac``), as the ENCODE ATAC-seq pipeline does it:

1. Reads: mapped, primary, not QC-failed, not duplicates (``-F 3844``), optionally ``MAPQ >= --mapq``. When the BAM has
   no duplicate marks (no ``samtools markdup`` / Picard ``MarkDuplicates`` in its header), duplicates are marked
   first (``samtools collate | fixmate -m | sort | markdup``; paired- or single-end).
2. Cut sites: ATAC-seq, the Tn5 insertion (5' end shifted +4 on the plus strand, -5 on the minus strand); DNase-seq,
   the 5' end. Each mate of a pair is a cut site.
3. ``macs2 callpeak -f BED --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01`` without a control,
   then ``macs2 bdgcmp -m ppois -S <cut sites / 1e6>``. The result is clipped to the BAM's chromosomes and sorted
   (bedGraph; also a bigWig when the output ends in .bw/.bigWig).

On osteoblast ATAC-seq the track agrees with the ENCODE ATAC pipeline's own p-value track (Spearman 0.996 over the
candidate elements). Applying the DNase procedure to ATAC-seq gave the same results as the ATAC procedure in osteoblasts,
so each assay keeps its own procedure. Needs ``samtools`` and ``bedtools`` (``pipeline`` pixi environment) and MACS2 2.2.7.1 (``macs2``
pixi environment; pass its path with ``--macs2``).
"""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile

ASSAYS = {"atac": (4, -5), "dnase": (0, 0)}   # 5'-end shift on the plus / minus strand
ROADMAP_FLAGS_EXCLUDED = 1796   # unmapped, secondary, QC fail, duplicate (Roadmap DNase tags)
ROADMAP_MAPQ, ROADMAP_TAG, ROADMAP_MAX_TAGS = 30, 36, 50_000_000
ROADMAP_CALLPEAK = "-p 0.01 --nomodel --extsize 146 -B --SPMR"
ROADMAP_CONTROL_URL = "https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/E123-DNaseControl.tagAlign.gz"
HG19_TO_HG38_CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
READ_FLAGS_EXCLUDED = 3844   # unmapped, secondary, QC fail, duplicate, supplementary
MACS2_CALLPEAK = "--nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01"


def _run(cmd: str, **kw) -> str:
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True, **kw).stdout


def has_duplicate_marks(bam, samtools: str = "samtools") -> bool:
    """True when the BAM header records a duplicate-marking step."""
    head = _run(f"{shlex.quote(samtools)} view -H {shlex.quote(str(bam))}")
    return any(line.startswith("@PG") and any(k in line for k in ("markdup", "MarkDuplicates", "rmdup")) for line in head.splitlines())


def is_paired(bam, samtools: str = "samtools", n: int = 10000) -> bool:
    flags = _run(f"set +o pipefail; {shlex.quote(samtools)} view {shlex.quote(str(bam))} | head -n {n} | cut -f2")   # head closes the pipe early
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


def roadmap_tags(bam, out, sizes: dict, samtools: str = "samtools", bedtools: str = "bedtools", threads: int = 4,
                 max_tags: int = ROADMAP_MAX_TAGS) -> int:
    """Roadmap DNase tags of a BAM (36-bp 5' tags, MAPQ >= 30, read 1 of pairs, <= ``max_tags``) as a 6-column BED."""
    f = f"-F {ROADMAP_FLAGS_EXCLUDED}" + (" -f 64" if is_paired(bam, samtools) else "")
    tmp = str(out) + ".all"
    with open(tmp + ".sizes", "w") as fh:
        fh.writelines(f"{c}\t{n}\n" for c, n in sizes.items())
    _run(f"{shlex.quote(samtools)} view -@ {int(threads)} -b {f} -q {ROADMAP_MAPQ} {shlex.quote(str(bam))}"
         f" | {shlex.quote(bedtools)} bamtobed -i stdin"
         f" | awk -v OFS='\\t' -v L={ROADMAP_TAG}"
         f" '{{if ($6 == \"+\") $3 = $2 + L; else {{$2 = $3 - L; if ($2 < 0) $2 = 0}}; print $1, $2, $3, \"N\", 1000, $6}}'"
         f" | awk 'NR == FNR {{L[$1] = $2; next}} ($1 in L) && $3 <= L[$1]' {shlex.quote(tmp + '.sizes')} - > {shlex.quote(tmp)}")
    n = int(_run(f"wc -l < {shlex.quote(tmp)}").strip())
    if n > max_tags:   # as in the reconstruction: shuf with a fixed random source
        _run(f"shuf -n {int(max_tags)} --random-source=<(yes) {shlex.quote(tmp)} > {shlex.quote(str(out))}")
        os.remove(tmp); n = max_tags
    else:
        os.replace(tmp, out)
    os.remove(tmp + ".sizes")
    return n


def roadmap_control(output, crossmap: str = "CrossMap", tmpdir=None, log=print) -> str:
    """Download Roadmap's uniform DNase control (E123-DNaseControl, hg19) and lift it to GRCh38 with CrossMap.
    Writes a 6-column BED (about 1.6 GB) for ``pvalue_track(..., assay='dnase', control=...)``."""
    import urllib.request
    if not shutil.which(crossmap):
        raise RuntimeError(f"{crossmap} not found (pixi install -e macs2)")
    tmp = tempfile.mkdtemp(dir=tmpdir or os.path.dirname(os.path.abspath(str(output))))
    try:
        ctrl, chain = os.path.join(tmp, "control.hg19.tagAlign.gz"), os.path.join(tmp, "hg19ToHg38.over.chain.gz")
        log(f"downloading {ROADMAP_CONTROL_URL}")
        urllib.request.urlretrieve(ROADMAP_CONTROL_URL, ctrl)
        urllib.request.urlretrieve(HG19_TO_HG38_CHAIN, chain)
        _run(f"zcat {shlex.quote(ctrl)} > {tmp}/control.hg19.bed")
        log("lifting to GRCh38 (CrossMap)")
        _run(f"{shlex.quote(crossmap)} bed {shlex.quote(chain)} {tmp}/control.hg19.bed {tmp}/control.hg38.bed")
        _run(f"awk -v OFS='\\t' '$1 !~ /_/ && $1 != \"chrM\" && $3 > $2' {tmp}/control.hg38.bed > {shlex.quote(str(output))}")
        log(f"control -> {output} ({int(_run(f'wc -l < {shlex.quote(str(output))}').strip())} tags)")
        return str(output)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def pvalue_track(bams, output, assay: str = "atac", mapq: int = 0, mark_dups: str = "auto", macs2: str = "macs2",
                 samtools: str = "samtools", bedtools: str = "bedtools", threads: int = 4, tmpdir=None, control=None,
                 log=print) -> dict:
    """Write the accessibility -log10 P track (bedGraph; also a bigWig when ``output`` ends in .bw/.bigWig) of one BAM
    file or of several (biological replicates, pooled: their tags are combined before the cap and MACS2, as Roadmap
    and ENCODE pool replicates). DNase-seq uses the Roadmap procedure and needs ``control`` (``roadmap_control``); ATAC-seq the ENCODE ATAC
    procedure (``mapq`` and ``mark_dups`` apply to it: 'auto' marks when the header shows no marking step, 'yes' or
    'no'). Returns a summary."""
    bams = [str(b) for b in ([bams] if isinstance(bams, (str, os.PathLike)) else bams)]
    assay = assay.lower()
    if assay not in ASSAYS:
        raise ValueError(f"assay must be one of {list(ASSAYS)}")
    if assay == "dnase" and not control:
        raise ValueError("DNase-seq needs Roadmap's uniform control: cislink roadmap-control -o control.bed, then --control control.bed")
    for tool in (samtools, bedtools, macs2):
        if not shutil.which(tool):
            raise RuntimeError(f"{tool} not found (MACS2: pixi install -e macs2, then --macs2 .pixi/envs/macs2/bin/macs2)")
    out = str(output); bigwig = out.lower().endswith((".bw", ".bigwig"))
    bedgraph = (out.rsplit(".", 1)[0] + ".bedGraph") if bigwig else out
    tmp = tempfile.mkdtemp(dir=tmpdir or os.path.dirname(os.path.abspath(out)))
    try:
        sizes = chrom_sizes(bams[0], samtools)
        if assay == "dnase":
            tags = os.path.join(tmp, "tags.bed")
            for i, b in enumerate(bams):   # each replicate's tags uncapped, then one cap over the pool
                k = roadmap_tags(b, os.path.join(tmp, f"tags{i}.bed"), sizes, samtools, bedtools, threads, max_tags=10**12)
                log(f"{b}: {k} DNase tags")
            n = _pool_and_cap([os.path.join(tmp, f"tags{i}.bed") for i in range(len(bams))], tags, ROADMAP_MAX_TAGS)
            if n == 0:
                raise RuntimeError("no reads pass the filters")
            nc = int(_run(f"wc -l < {shlex.quote(str(control))}").strip())
            log(f"{n} DNase tags (Roadmap procedure), control {nc}; MACS2")
            m, t = shlex.quote(macs2), shlex.quote(tmp)
            _run(f"{m} callpeak -t {shlex.quote(tags)} -c {shlex.quote(str(control))} -f BED -n acc -g hs --outdir {t} {ROADMAP_CALLPEAK}")
            _run(f"{m} bdgcmp -t {t}/acc_treat_pileup.bdg -c {t}/acc_control_lambda.bdg -m ppois -S {min(n, nc) / 1e6!r} --o-prefix {t}/acc")
            return _finish(tmp, sizes, bedgraph, out, bigwig, {"tags": n, "assay": assay}, log)
        tags = os.path.join(tmp, "cuts.bed")
        for i, b in enumerate(bams):   # each replicate deduplicated on its own, cut sites pooled
            src = b
            if mark_dups == "yes" or (mark_dups == "auto" and not has_duplicate_marks(b, samtools)):
                log(f"{b}: marking duplicates (samtools fixmate/markdup)")
                src = mark_duplicates(b, os.path.join(tmp, f"markdup{i}.bam"), samtools, threads, tmp)
            k = cut_sites(src, os.path.join(tmp, f"cuts{i}.bed"), assay, mapq, samtools, bedtools, threads)
            log(f"{b}: {k} cut sites")
        n = _pool_and_cap([os.path.join(tmp, f"cuts{i}.bed") for i in range(len(bams))], tags, None)
        if n == 0:
            raise RuntimeError("no reads pass the filters")
        log(f"{n} cut sites ({assay}); MACS2")
        m, t = shlex.quote(macs2), shlex.quote(tmp)
        _run(f"{m} callpeak -t {shlex.quote(tags)} -f BED -n acc -g hs --outdir {t} {MACS2_CALLPEAK}")
        _run(f"{m} bdgcmp -t {t}/acc_treat_pileup.bdg -c {t}/acc_control_lambda.bdg -m ppois -S {n / 1e6!r} --o-prefix {t}/acc")
        return _finish(tmp, sizes, bedgraph, out, bigwig, {"cut_sites": n, "assay": assay}, log)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _pool_and_cap(parts, out, max_tags) -> int:
    """Concatenate tag files (removing them) and keep at most ``max_tags`` (fixed-seed shuf); returns the tag count."""
    _run(f"cat {' '.join(shlex.quote(p) for p in parts)} > {shlex.quote(str(out))}")
    for p in parts:
        os.remove(p)
    n = int(_run(f"wc -l < {shlex.quote(str(out))}").strip())
    if max_tags and n > max_tags:
        q = shlex.quote(str(out))
        _run(f"shuf -n {int(max_tags)} --random-source=<(yes) {q} > {q}.cap && mv {q}.cap {q}")
        n = max_tags
    return n


def _finish(tmp, sizes, bedgraph, out, bigwig, summary, log) -> dict:
    """Clip the MACS2 -log10 P bedGraph to the BAM's chromosomes, sort it and write the bigWig when asked."""
    t = shlex.quote(tmp)
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
    return {**summary, "bedgraph": bedgraph, "bigwig": out if bigwig else None}


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
