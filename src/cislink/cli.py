"""Command line: ``cislink pvalue``, ``cislink roadmap-control``, ``cislink reference download``, ``cislink validate``,
``cislink build``, ``cislink add-f1``, ``cislink score``, ``cislink rank`` and ``cislink info``."""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from . import __version__
from .model import MODELS, load_model


def _score(a) -> int:
    t = pd.read_csv(a.input, sep=a.sep, low_memory=False)
    for m in a.model:
        model = load_model(m)
        t[model.name] = model.predict(t)
        if a.threshold:
            t[f"{model.name} >= {a.threshold}% recall threshold"] = t[model.name] >= model.threshold(a.threshold)
    keep = [c for c in a.keep.split(",") if c] if a.keep else None
    if keep:
        t = t[keep + [c for c in t.columns if c.startswith("CisLink-")]]
    t.to_csv(a.output or sys.stdout, sep="\t", index=False, float_format="%.6g")
    return 0


def _load_inputs(a):
    import yaml

    from .inputs import read_variants
    cfg = yaml.safe_load(open(a.config))
    v = read_variants(a.variants) if a.variants else None
    e = None
    if a.elements:
        e = pd.read_csv(a.elements, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3], names=["chr", "start", "end", "name"])
    return cfg, v, e


def _validate(a) -> int:
    from .inputs import InputError, preflight
    cfg, v, e = _load_inputs(a)
    try:
        preflight(cfg, v, e, log=lambda s: print(s, file=sys.stderr))
    except InputError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    return 0


def _build(a) -> int:
    from . import features as cf
    from .inputs import InputError, preflight
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
    from . import roadmap as rm
    cfg, v, e = _load_inputs(a)
    try:
        preflight(cfg, v, e, log=log)
        ref = rm.reference_path(a.roadmap_ref or cfg.get("roadmap_reference"))
    except (InputError, rm.ReferenceError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    crossmap = a.crossmap or cfg.get("crossmap", "CrossMap")
    tx, _ = cf.read_gtf(cfg["gtf"])
    tss = cf.protein_coding_tss(tx)
    rows = cf.variant_pairs(v, tss, a.max_distance) if v is not None else cf.element_pairs(e, tss, a.max_distance)
    bt, th = cfg.get("bedtools", "bedtools"), int(cfg.get("threads", 4))
    # the cell's quantiles: read when given, otherwise computed and saved next to the output for reuse
    if cfg.get("chromatin_quantiles"):
        cq = pd.read_csv(cfg["chromatin_quantiles"], sep="\t")
    else:
        candidates = cfg.get("candidates") or cf.make_candidates(cfg["accessibility_pvalue"], a.output + ".candidates.bed", bt)
        log("cell chromatin quantiles over candidate elements")
        cq = cf.cell_chromatin_quantiles(candidates, cfg["tracks"], bt, th)
        cq.to_csv(a.output + ".chromatin_quantiles.tsv", sep="\t", index=False, float_format="%.6g")
    hq = None
    map_hic = cfg.get("map_hic")
    if map_hic is None:   # default: map only maps finer than the training 5 kb
        map_hic = int(cfg.get("hic_resolution", cf.HIC_BIN)) < cf.HIC_BIN
    if cfg.get("hic") is not None and bool(map_hic):
        if cfg.get("hic_quantiles"):
            hq = pd.read_csv(cfg["hic_quantiles"], sep="\t")
        else:
            log("cell Hi-C quantiles (20 random 2-Mb windows per autosome)")
            hq = cf.hic_cell_quantiles(cfg["hic"], int(cfg.get("hic_resolution", cf.HIC_BIN)), cfg.get("hic_norm", "KR"))
            hq.to_csv(a.output + ".hic_quantiles.tsv", sep="\t", index=False, float_format="%.6g")
    out = cf.build_features(rows, gtf=cfg["gtf"], tracks=cfg["tracks"], transcript_quantification=cfg["transcript_quantification"],
                            gene_quantification=cfg["gene_quantification"], candidates_bed=None, cell_quantiles=cq, hic_quantiles=hq, hic=cfg.get("hic"),
                            hic_resolution=int(cfg.get("hic_resolution", cf.HIC_BIN)), hic_norm=cfg.get("hic_norm", "KR"),
                            map_hic=bool(map_hic), megamap=cfg.get("megamap", cf.MEGAMAP),
                            bedtools=bt, threads=th, log=log)
    out = rm.add_f1(out, cfg["gtf"], ref, crossmap=crossmap, log=log)
    if "SNPID" in rows:
        out.insert(1, "SNPID", rows.SNPID.to_numpy())
    else:
        out.insert(1, "element", rows.element.to_numpy())
    if a.score:
        for m in MODELS:
            model = load_model(m)
            if m == "e2g" and cfg.get("hic") is None:
                continue
            out[model.name] = model.predict(out)
    out.to_csv(a.output, sep="\t", index=False, float_format="%.8g")
    print(f"wrote {a.output}: {len(out)} rows", file=sys.stderr)
    return 0


def _add_f1(a) -> int:
    from . import roadmap as rm
    try:
        ref = rm.reference_path(a.roadmap_ref)
    except rm.ReferenceError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    t = pd.read_csv(a.input, sep="\t", low_memory=False)
    need = {"Enh", "TX"} - set(t.columns)
    if need:
        print(f"error: the table needs columns Enh (chr:start:end) and TX (gene optional); missing {sorted(need)}", file=sys.stderr)
        return 1
    t = rm.add_f1(t, a.gtf, ref, crossmap=a.crossmap, genome=a.genome, log=lambda s: print(s, file=sys.stderr, flush=True))
    t.to_csv(a.output, sep="\t", index=False, float_format="%.8g")
    return 0


def _reference(a) -> int:
    from . import roadmap as rm
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
    try:
        if a.action == "download":
            ref = rm.download(a.dest, base_url=a.base_url, keep_archives=a.keep_archives, log=log)
            print(f"reference ready: {ref}\nuse it with: cislink build --roadmap-ref {ref}   (or export {rm.ENV}={ref})")
        else:
            log(f"verified {rm.verify(rm.reference_path(a.dest))} files")
    except rm.ReferenceError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    return 0


def _pvalue(a) -> int:
    from .pvalue import pvalue_track
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
    pvalue_track(a.bam, a.output, assay=a.assay, mapq=a.mapq, mark_dups=a.mark_duplicates, macs2=a.macs2,
                 samtools=a.samtools, bedtools=a.bedtools, threads=a.threads, control=a.control, log=log)
    return 0


def _roadmap_control(a) -> int:
    from .pvalue import roadmap_control
    roadmap_control(a.output, crossmap=a.crossmap, log=lambda s: print(s, file=sys.stderr, flush=True))
    return 0


def _rank(a) -> int:
    """Pair scores (maximum over a gene's transcript rows) and gene ranks per variant or element."""
    t = pd.read_csv(a.input, sep="\t", low_memory=False)
    key = "SNPID" if "SNPID" in t else "element"
    cols = [c for c in ("CisLink-E2G", "CisLink-V2G") if c in t]
    if not cols:
        raise SystemExit("no CisLink score columns; run 'cislink score' or 'cislink build --score' first")
    p = t.groupby([key, "gene"], sort=False).agg({**{c: "max" for c in cols}, "distance": "min"}).reset_index()
    for c in cols:
        model = load_model(c)
        p[f"{c} rank"] = p.groupby(key)[c].rank(ascending=False, method="min").astype(int)
        p[f"{c} >= 70% recall threshold"] = p[c] >= model.threshold(70)
    p = p.sort_values([key, f"{cols[0]} rank"])
    try:
        p.to_csv(a.output or sys.stdout, sep="\t", index=False, float_format="%.6g")
    except BrokenPipeError:   # e.g. piped into head
        pass
    return 0


def _info(a) -> int:
    for m in MODELS:
        model = load_model(m)
        print(f"{model.name}: {len(model.features)} features; thresholds {model.thresholds.get('recall_70')} (70% recall), "
              f"{model.thresholds.get('recall_50')} (50% recall)")
        print("  " + ", ".join(model.features))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="cislink", description="Score element-gene and variant-gene pairs with CisLink.")
    p.add_argument("--version", action="version", version=f"cislink {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score", help="score a harmonized feature table")
    s.add_argument("input", help="feature table (TSV by default) with the model's feature columns")
    s.add_argument("-m", "--model", nargs="+", default=list(MODELS), choices=MODELS)
    s.add_argument("-o", "--output", help="output TSV (default stdout)")
    s.add_argument("--sep", default="\t")
    s.add_argument("--keep", help="comma-separated columns to keep with the scores (default: all)")
    s.add_argument("--threshold", type=int, choices=(50, 70), help="add a column flagging scores above the recall threshold")
    s.set_defaults(func=_score)
    def inputs_args(b):
        g = b.add_mutually_exclusive_group(required=True)
        g.add_argument("--variants", help="TSV with columns SNPID, chr, pos (1-based, GRCh38), or chrom, pos[, ID] without a header")
        g.add_argument("--elements", help="BED (chrom, start, end, name) of candidate elements (GRCh38)")
        b.add_argument("--config", required=True, help="YAML with the cell's inputs")
    vd = sub.add_parser("validate", help="check a build configuration and its inputs without running the pipeline")
    inputs_args(vd); vd.set_defaults(func=_validate)
    b = sub.add_parser("build", help="compute harmonized features from a cell type's raw data (see docs/pipeline.md)")
    inputs_args(b)
    b.add_argument("--max-distance", type=int, default=1_000_000)
    b.add_argument("--score", action="store_true", help="also add CisLink-V2G (and CisLink-E2G when the config has Hi-C) scores")
    b.add_argument("--roadmap-ref", help="Roadmap activity reference directory (cislink reference download); default $CISLINK_ROADMAP_REF "
                                         "or 'roadmap_reference' in the config")
    b.add_argument("--crossmap", help="CrossMap executable for the GRCh38 -> hg19 lift of the F1 features (default: config 'crossmap' or CrossMap)")
    b.add_argument("-o", "--output", required=True)
    b.set_defaults(func=_build)
    f1 = sub.add_parser("add-f1", help="add the F1 features (Roadmap activity-expression correlations) to a feature table")
    f1.add_argument("input", help="feature table with columns Enh (chr:start:end), TX and optionally gene (symbol)")
    f1.add_argument("--gtf", required=True, help="the annotation used to build the table (transcript -> gene IDs)")
    f1.add_argument("-o", "--output", required=True)
    f1.add_argument("--roadmap-ref", help="reference directory (default $CISLINK_ROADMAP_REF)")
    f1.add_argument("--crossmap", default="CrossMap"); f1.add_argument("--genome", choices=("GRCh38", "hg19"), default="GRCh38")
    f1.set_defaults(func=_add_f1)
    rf = sub.add_parser("reference", help="download or verify the Roadmap activity reference (16.3 GB) used by the F1 features")
    rf.add_argument("action", choices=("download", "verify"))
    rf.add_argument("--dest", required=True, help="download: directory to put it in; verify: the reference directory")
    rf.add_argument("--base-url", help="location of the archives (default: the Zenodo record)")
    rf.add_argument("--keep-archives", action="store_true")
    rf.set_defaults(func=_reference)
    pv = sub.add_parser("pvalue", help="accessibility -log10 P track from an ATAC-seq or DNase-seq BAM (see docs/pipeline.md)")
    pv.add_argument("bam", nargs="+", help="BAM file(s); biological replicates are pooled")
    pv.add_argument("-o", "--output", required=True, help="bedGraph, or .bw/.bigWig (needs pyBigWig)")
    pv.add_argument("--assay", choices=("atac", "dnase"), required=True,
                    help="atac: ENCODE ATAC procedure (Tn5 insertions); dnase: Roadmap procedure (36-bp tags, uniform control)")
    pv.add_argument("--control", help="DNase only: Roadmap uniform control BED on GRCh38 (cislink roadmap-control)")
    pv.add_argument("--mapq", type=int, default=0, help="ATAC only: minimum mapping quality (default 0: the BAM's own filtering)")
    pv.add_argument("--mark-duplicates", choices=("auto", "yes", "no"), default="auto",
                    help="ATAC only: mark duplicates first: auto = when the BAM header shows no marking step")
    pv.add_argument("--macs2", default="macs2", help="MACS2 2.2.7.1 executable (pixi install -e macs2: .pixi/envs/macs2/bin/macs2)")
    pv.add_argument("--threads", type=int, default=4); pv.add_argument("--samtools", default="samtools"); pv.add_argument("--bedtools", default="bedtools")
    pv.set_defaults(func=_pvalue)
    rc = sub.add_parser("roadmap-control", help="download Roadmap's uniform DNase control and lift it to GRCh38 (for pvalue --assay dnase)")
    rc.add_argument("-o", "--output", required=True, help="output BED (about 1.6 GB)")
    rc.add_argument("--crossmap", default="CrossMap", help="CrossMap executable (pixi install -e macs2: .pixi/envs/macs2/bin/CrossMap)")
    rc.set_defaults(func=_roadmap_control)
    r = sub.add_parser("rank", help="pair scores and gene ranks per variant or element from a scored table")
    r.add_argument("input"); r.add_argument("-o", "--output", help="output TSV (default stdout)")
    r.set_defaults(func=_rank)
    i = sub.add_parser("info", help="list the models' features and thresholds")
    i.set_defaults(func=_info)
    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
