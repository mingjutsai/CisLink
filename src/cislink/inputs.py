"""Input handling for ``cislink build``: chromosome naming, expression tables and a preflight check.

The preflight turns input problems into clear errors before any long step runs. The most important check is
chromosome naming: a bedGraph that writes ``1`` while the regions say ``chr1`` gives all-zero signal without any
error, so naming is detected per file and converted (and tracks that disagree with each other are rejected).
"""
from __future__ import annotations

import gzip
import os
import shutil

import pandas as pd

MARKS = ("DNase", "H4K20me1", "H3K9me3", "H3K9ac", "H3K79me2", "H3K4me3", "H3K4me2", "H3K4me1", "H3K36me3",
         "H3K27me3", "H3K27ac", "H2A.Z", "CTCF")


def normalize_ensembl_id(identifier: str) -> str:
    """Stable Ensembl ID: ``ENSG00000123456.11`` -> ``ENSG00000123456``."""
    return str(identifier).strip().split(".", 1)[0]


def strip_chr(c: str) -> str:
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def format_chr(c: str, use_prefix: bool) -> str:
    return ("chr" + strip_chr(c)) if use_prefix else strip_chr(c)


def is_bigwig(path) -> bool:
    return str(path).lower().endswith((".bw", ".bigwig"))


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _rows(path, n: int = 1):
    """The first ``n`` data rows (skipping blank, '#', 'track' and 'browser' lines)."""
    out = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            out.append(line.rstrip("\n").split("\t"))
            if len(out) >= n:
                break
    return out


def chr_prefix(path) -> bool | None:
    """True when a file's chromosomes are named ``chrN``, False for ``N``; None for an empty file. bigWig and .hic
    files are read through their own libraries."""
    if is_bigwig(path):
        import pyBigWig
        bw = pyBigWig.open(str(path)); names = list(bw.chroms()); bw.close()
        return any(n.lower().startswith("chr") for n in names) if names else None
    rows = _rows(path)
    return rows[0][0].lower().startswith("chr") if rows else None


def read_expression_raw(path, kind: str) -> pd.DataFrame:
    """Every row of an expression table as ``id`` (stable Ensembl ID) and ``TPM``, without filtering (for scale factors)."""
    idcol = "transcript_id" if kind == "transcript" else "gene_id"
    head = pd.read_csv(path, sep="\t", nrows=0)
    if idcol in head.columns and "TPM" in head.columns:
        t = pd.read_csv(path, sep="\t", usecols=[idcol, "TPM"]).rename(columns={idcol: "id"})
    else:
        t = pd.read_csv(path, sep="\t", header=None, comment="#").iloc[:, [0, -1]].set_axis(["id", "TPM"], axis=1)
        t = t[pd.to_numeric(t.TPM, errors="coerce").notna()]
    return pd.DataFrame({"id": t.id.map(normalize_ensembl_id), "TPM": t.TPM.astype(float)})


def read_expression(path, kind: str) -> pd.DataFrame:
    """Expression table as columns ``id`` (stable Ensembl ID) and ``TPM``. ``kind``: 'transcript' or 'gene'.

    Accepted: ENCODE RSEM output (``transcript_id`` / ``gene_id`` and ``TPM`` columns), or any table with the
    Ensembl ID in the first column and the value in the last (with or without a header). PAR_Y entries are dropped;
    for duplicated stable IDs the last row wins."""
    idcol = "transcript_id" if kind == "transcript" else "gene_id"
    head = pd.read_csv(path, sep="\t", nrows=0, comment=None)
    if idcol in head.columns and "TPM" in head.columns:
        t = pd.read_csv(path, sep="\t", usecols=[idcol, "TPM"]).rename(columns={idcol: "id"})
    else:
        t = pd.read_csv(path, sep="\t", header=None, comment="#")
        t = t.iloc[:, [0, -1]].set_axis(["id", "TPM"], axis=1)
        t = t[pd.to_numeric(t.TPM, errors="coerce").notna()]   # drop a header line, if any
    t = t[~t.id.astype(str).str.contains("_PAR_Y")].copy()
    t["id"] = t.id.map(normalize_ensembl_id); t["TPM"] = t.TPM.astype(float)
    return t.drop_duplicates("id", keep="last").reset_index(drop=True)


def read_variants(path) -> pd.DataFrame:
    """Variants as SNPID, chr, pos. A header with these names is used; otherwise columns are chromosome, 1-based
    position and (optionally) ID, as in a BED-like table."""
    head = pd.read_csv(path, sep="\t", nrows=0)
    if {"SNPID", "chr", "pos"} <= set(head.columns):
        v = pd.read_csv(path, sep="\t", dtype={"SNPID": str, "chr": str})
    else:
        raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
        v = pd.DataFrame({"chr": raw[0], "pos": raw[1]})
        v["SNPID"] = raw[2] if raw.shape[1] > 2 else v.chr.map(strip_chr) + ":" + v.pos
    v["pos"] = pd.to_numeric(v.pos, errors="raise").astype(int)
    return v[["SNPID", "chr", "pos"]]


class InputError(ValueError):
    pass


def preflight(cfg: dict, variants: pd.DataFrame | None = None, elements: pd.DataFrame | None = None, log=print) -> dict:
    """Check a ``cislink build`` configuration and its inputs; return the detected layout. Raises InputError."""
    need = ["gtf", "transcript_quantification", "gene_quantification", "tracks"]
    missing = [k for k in need if not cfg.get(k)]
    if missing:
        raise InputError("missing config keys: " + ", ".join(missing))
    for k in ("gtf", "transcript_quantification", "gene_quantification"):
        if not os.path.isfile(cfg[k]):
            raise InputError(f"{k} is not a readable file: {cfg[k]}")
    if not (cfg.get("chromatin_quantiles") or cfg.get("candidates") or cfg.get("accessibility_pvalue")):
        raise InputError("give 'candidates' (BED), 'accessibility_pvalue' (bedGraph to call them) or 'chromatin_quantiles'")

    tracks = cfg["tracks"]
    unknown = [m for m in tracks if m not in MARKS]
    if unknown:
        raise InputError(f"unknown track names {unknown}; use {list(MARKS)}")
    if "DNase" not in tracks:
        raise InputError("an accessibility track ('DNase': ATAC-seq or DNase-seq) is required")
    prefixes = {}
    for m, p in tracks.items():
        if not os.path.isfile(p):
            raise InputError(f"track {m} is not a readable file: {p}")
        if not is_bigwig(p):
            row = _rows(p)
            if not row or len(row[0]) < 4:
                raise InputError(f"track {m}: bedGraph rows need chrom, start, end, value ({p})")
            try:
                int(row[0][1]); int(row[0][2]); float(row[0][3])
            except ValueError:
                raise InputError(f"track {m}: first data row is not chrom, start, end, value: {row[0][:4]} ({p})") from None
        prefixes[m] = chr_prefix(p)
    bedgraph_prefix = {prefixes[m] for m, p in tracks.items() if not is_bigwig(p)}
    if len(bedgraph_prefix) > 1:
        bad = {m: ("chrN" if v else "N") for m, v in prefixes.items() if not is_bigwig(tracks[m])}
        raise InputError(f"bedGraph tracks disagree on chromosome naming: {bad}")
    bt = cfg.get("bedtools", "bedtools")
    if (bedgraph_prefix or cfg.get("accessibility_pvalue")) and not shutil.which(bt):
        raise InputError(f"bedtools not found ({bt}); it is needed for bedGraph tracks")

    # identifiers shared by the annotation and the expression tables
    from .features import read_gtf
    tx, genes = read_gtf(cfg["gtf"])
    gtf_prefix = tx.chrom.iloc[0].lower().startswith("chr") if len(tx) else True
    tq, gq = read_expression(cfg["transcript_quantification"], "transcript"), read_expression(cfg["gene_quantification"], "gene")
    tx_hit = set(tx.TX.map(normalize_ensembl_id)) & set(tq.id)
    g_hit = set(genes.gene_id.map(normalize_ensembl_id)) & set(gq.id)
    if not tx_hit or not g_hit:
        raise InputError("the annotation and expression files share no stable Ensembl IDs "
                         f"(transcripts {len(tx_hit)}, genes {len(g_hit)})")
    pc = tx[tx.coding]
    miss = 1 - pc.TX.map(normalize_ensembl_id).isin(tq.id).mean()
    if miss > 0.2:
        log(f"warning: {miss:.0%} of protein-coding transcripts have no expression entry (scored as 0)")

    if variants is not None:
        if variants.SNPID.duplicated().any():
            raise InputError(f"duplicated variant IDs, e.g. {variants.SNPID[variants.SNPID.duplicated()].iloc[0]}")
        if (variants.pos < 1).any():
            raise InputError("variant positions must be 1-based and positive")
        unknown = set(variants.chr.map(lambda c: format_chr(c, gtf_prefix))) - set(tx.chrom)
        if unknown:
            log(f"warning: chromosomes not in the annotation (no candidate genes): {sorted(unknown)[:5]}")
    if elements is not None and (elements.end <= elements.start).any():
        raise InputError("elements need end > start")
    for k in ("hic", "megamap"):
        p = cfg.get(k)
        if p and not str(p).startswith(("http://", "https://")) and not os.path.isfile(p):
            raise InputError(f"{k} is not a readable file: {p}")

    style = lambda v: "chrN" if v else "N"  # noqa: E731
    info = {"gtf_chr_prefix": gtf_prefix, "track_chr_prefix": prefixes, "transcripts_matched": len(tx_hit), "genes_matched": len(g_hit)}
    log(f"preflight passed: {len(tracks)} tracks ({', '.join(f'{m} {style(v)}' for m, v in prefixes.items())}); "
        f"annotation {style(gtf_prefix)}; {len(tx_hit)} transcripts and {len(g_hit)} genes with expression")
    return info
