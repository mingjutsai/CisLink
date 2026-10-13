"""Command line: ``cislink build``, ``cislink score``, ``cislink rank`` and ``cislink info``."""
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


def _build(a) -> int:
    import yaml

    from . import features as cf
    cfg = yaml.safe_load(open(a.config))
    tx, _ = cf.read_gtf(cfg["gtf"])
    tss = cf.protein_coding_tss(tx)
    if a.variants:
        v = pd.read_csv(a.variants, sep="\t", dtype={"SNPID": str})
        rows = cf.variant_pairs(v, tss, a.max_distance)
    else:
        e = pd.read_csv(a.elements, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "name"])
        rows = cf.element_pairs(e, tss, a.max_distance)
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
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
    if cfg.get("hic") is not None and bool(cfg.get("map_hic", True)):
        if cfg.get("hic_quantiles"):
            hq = pd.read_csv(cfg["hic_quantiles"], sep="\t")
        else:
            log("cell Hi-C quantiles (20 random 2-Mb windows per autosome)")
            hq = cf.hic_cell_quantiles(cfg["hic"], int(cfg.get("hic_resolution", cf.HIC_BIN)), cfg.get("hic_norm", "KR"))
            hq.to_csv(a.output + ".hic_quantiles.tsv", sep="\t", index=False, float_format="%.6g")
    out = cf.build_features(rows, gtf=cfg["gtf"], tracks=cfg["tracks"], transcript_quantification=cfg["transcript_quantification"],
                            gene_quantification=cfg["gene_quantification"], candidates_bed=None, cell_quantiles=cq, hic_quantiles=hq, hic=cfg.get("hic"),
                            hic_resolution=int(cfg.get("hic_resolution", cf.HIC_BIN)), hic_norm=cfg.get("hic_norm", "KR"),
                            map_hic=bool(cfg.get("map_hic", True)), megamap=cfg.get("megamap", cf.MEGAMAP),
                            bedtools=bt, threads=th, log=log)
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
    b = sub.add_parser("build", help="compute harmonized features from a cell type's raw data (see docs/pipeline.md)")
    g = b.add_mutually_exclusive_group(required=True)
    g.add_argument("--variants", help="TSV with columns SNPID, chr (chrN), pos (1-based, GRCh38)")
    g.add_argument("--elements", help="BED (chrom, start, end, name) of candidate elements (GRCh38)")
    b.add_argument("--config", required=True, help="YAML with the cell's inputs")
    b.add_argument("--max-distance", type=int, default=1_000_000)
    b.add_argument("--score", action="store_true", help="also add CisLink-V2G (and CisLink-E2G when the config has Hi-C) scores")
    b.add_argument("-o", "--output", required=True)
    b.set_defaults(func=_build)
    r = sub.add_parser("rank", help="pair scores and gene ranks per variant or element from a scored table")
    r.add_argument("input"); r.add_argument("-o", "--output", help="output TSV (default stdout)")
    r.set_defaults(func=_rank)
    i = sub.add_parser("info", help="list the models' features and thresholds")
    i.set_defaults(func=_info)
    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
