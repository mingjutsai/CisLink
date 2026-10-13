"""Compute CisLink input features from a cell type's raw data (GRCh38).

The steps reproduce the pipeline used for the manuscript:

1. **Pairs.** Every variant (or element) against every protein-coding transcript whose TSS lies within 1 Mb.
   A variant's element is the variant +/-250 bp; the promoter is the TSS +/-250 bp. One row per transcript;
   a pair's score is the maximum over its rows.
2. **Raw signals.** For each chromatin track, sum(signal x overlapping bp) over the element and the promoter
   (bedGraph through ``bedtools intersect``, or bigWig through pyBigWig). Expression: transcript and gene TPM.
3. **Harmonization to K562** (``cislink.harmonize``): chromatin quantile-mapped over the cell's genome-wide candidate
   elements; expression as protein-coding TPM; Hi-C observed/expected mapped within distance strata.
4. **Distance and gene context**: distance, log distance, element length, power-law contact, genes between element
   and TSS, distance rank of the gene, nearest-gene status, ubiquitously expressed gene.
5. **Hi-C**: the cell's own map (CisLink-E2G, ``HiC_Contact``) and the multi-tissue megamap (CisLink-V2G,
   ``HiC_megamap_oe``), observed/expected between the element-midpoint bin and the TSS bin.

Hi-C needs ``hic-straw``; bigWig tracks need ``pyBigWig``; bedGraph tracks and candidate calling need ``bedtools``.
"""
from __future__ import annotations

import gzip
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from importlib import resources

import numpy as np
import pandas as pd

from . import harmonize as hz

MARKS = hz.MARKS
ELEMENT_HALF_WIDTH = 250
PROMOTER_HALF_WIDTH = 250
HIC_BIN = 5000
POWERLAW_GAMMA = 1.024
MAX_DISTANCE = 1_000_000
MEGAMAP = "https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic"
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CANDIDATE_THRESHOLDS = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 8, 10]
CANDIDATE_TARGET = 150_000
CANDIDATE_MERGE_GAP = 100


def _sig6(q: pd.DataFrame) -> pd.DataFrame:
    """Quantiles at 6 significant digits, as the manuscript pipeline stored and re-read them. The top of a distribution
    is steep, so this matters for a few high-signal promoters; keeping it makes results reproducible exactly."""
    return q.apply(lambda c: c.map(lambda v: float(f"{v:.6g}") if np.isfinite(v) else v))


# ---------------------------------------------------------------- annotation
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gtf(gtf) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Transcripts (chrom, tss, TX, gene, coding, transcript_type) in file order, and genes (chrom, tss, gene, gene_id,
    gene_type). ``coding``: 'protein_coding' appears anywhere on the line (the rule for candidate TSSs)."""
    tx, genes = [], []
    with _open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("transcript", "gene"):
                continue
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            info = f[8]
            name = re.search(r'gene_name "([^"]+)"', info).group(1)
            if f[2] == "transcript":
                tid = re.sub(r"_\d+$", "", re.search(r'transcript_id "([^"]+)"', info).group(1))
                tt = re.search(r'transcript_type "([^"]+)"', info)
                tx.append((f[0], tss, tid, name, "protein_coding" in line, tt.group(1) if tt else ""))
            else:
                genes.append((f[0], tss, name, re.search(r'gene_id "([^"]+)"', info).group(1),
                              re.search(r'gene_type "([^"]+)"', info).group(1)))
    return (pd.DataFrame(tx, columns=["chrom", "tss", "TX", "gene", "coding", "transcript_type"]),
            pd.DataFrame(genes, columns=["chrom", "tss", "gene", "gene_id", "gene_type"]))


def protein_coding_tss(transcripts: pd.DataFrame) -> pd.DataFrame:
    """Protein-coding transcript TSSs; a transcript whose TSS equals that of the preceding one is dropped."""
    t = transcripts[transcripts.coding]
    return t[t.tss.ne(t.tss.shift())].drop(columns=["coding", "transcript_type"]).reset_index(drop=True)


def ubiquitous_genes() -> set[str]:
    """ABC's list of ubiquitously expressed genes (Fulco et al. 2019)."""
    with resources.files("cislink.reference").joinpath("UbiquitouslyExpressedGenes.txt").open() as fh:
        return set(fh.read().split())


# ---------------------------------------------------------------- pairs
def variant_pairs(variants: pd.DataFrame, tss: pd.DataFrame, max_distance: int = MAX_DISTANCE) -> pd.DataFrame:
    """Rows of variant x protein-coding transcript with the TSS within ``max_distance``.
    ``variants``: columns SNPID, chr (chrN), pos (1-based)."""
    rows = []
    for c, v in variants.drop_duplicates("SNPID").groupby("chr", sort=False):
        t = tss[tss.chrom == c].sort_values("tss", kind="mergesort")
        pos = t.tss.to_numpy()
        for snp, p in zip(v.SNPID, v.pos):
            lo, hi = np.searchsorted(pos, p - max_distance, "left"), np.searchsorted(pos, p + max_distance, "right")
            sub = t.iloc[lo:hi]
            rows.append(pd.DataFrame({"SNPID": snp, "chrom": c, "pos": p, "TX": sub.TX.to_numpy(), "gene": sub.gene.to_numpy(),
                                      "tss": sub.tss.to_numpy()}))
    r = pd.concat(rows, ignore_index=True)
    r["es"], r["ee"] = r.pos - ELEMENT_HALF_WIDTH, r.pos + ELEMENT_HALF_WIDTH
    r["ps"], r["pe"] = r.tss - PROMOTER_HALF_WIDTH, r.tss + PROMOTER_HALF_WIDTH
    r["ID"] = r.SNPID.astype(str) + ":" + r.gene + "|" + r.TX
    return r


def element_pairs(elements: pd.DataFrame, tss: pd.DataFrame, max_distance: int = MAX_DISTANCE) -> pd.DataFrame:
    """Rows of element x protein-coding transcript with the TSS within ``max_distance`` of the element midpoint.
    ``elements``: columns name, chr, start, end (0-based half-open, as BED)."""
    rows = []
    for c, v in elements.groupby("chr", sort=False):
        t = tss[tss.chrom == c].sort_values("tss", kind="mergesort")
        pos = t.tss.to_numpy()
        for name, s, e in zip(v.name, v.start, v.end):
            m = (s + e) // 2
            lo, hi = np.searchsorted(pos, m - max_distance, "left"), np.searchsorted(pos, m + max_distance, "right")
            sub = t.iloc[lo:hi]
            rows.append(pd.DataFrame({"element": name, "chrom": c, "es": s, "ee": e, "TX": sub.TX.to_numpy(),
                                      "gene": sub.gene.to_numpy(), "tss": sub.tss.to_numpy()}))
    r = pd.concat(rows, ignore_index=True)
    r["ps"], r["pe"] = r.tss - PROMOTER_HALF_WIDTH, r.tss + PROMOTER_HALF_WIDTH
    r["ID"] = r.element.astype(str) + ":" + r.gene + "|" + r.TX
    return r


# ---------------------------------------------------------------- chromatin signal
def _is_bigwig(path) -> bool:
    return str(path).lower().endswith((".bw", ".bigwig"))


def _sums_bigwig(regions: pd.DataFrame, path) -> np.ndarray:
    import pyBigWig
    bw = pyBigWig.open(str(path)); chroms = bw.chroms(); out = np.zeros(len(regions))
    for i, (c, s, e) in enumerate(zip(regions.chrom, regions.start, regions.end)):
        if c not in chroms:
            continue
        iv = bw.intervals(c, max(int(s), 0), min(int(e), chroms[c])) or []
        out[i] = sum(v * (min(e, b) - max(s, a)) for a, b, v in iv if min(e, b) - max(s, a) > 0)
    bw.close()
    return out


def _sums_bedgraph(n_regions: int, bed: str, path, bedtools: str, tmp: str, tag: str) -> np.ndarray:
    out = os.path.join(tmp, f"{tag}.sums")
    # bedtools reports an overlap of -1 bp for a zero-length interval on a boundary; count it as none
    cmd = (f"set -o pipefail; '{bedtools}' intersect -a '{path}' -b '{bed}' -wo"
           f" | awk -F'\\t' '$9 > 0 {{s[$8] += $4 * $9}} END {{for (k in s) printf \"%s\\t%.6f\\n\", k, s[k]}}' > '{out}'")
    subprocess.run(["bash", "-c", cmd], check=True)
    s = pd.read_csv(out, sep="\t", header=None, names=["i", "v"]).set_index("i").v
    return s.reindex(np.arange(n_regions)).fillna(0.0).to_numpy(dtype=float)


def signal_sums(regions: pd.DataFrame, tracks: dict[str, str], bedtools: str = "bedtools", threads: int = 4) -> pd.DataFrame:
    """sum(signal x overlapping bp) per region (columns chrom, start, end) for each track; NaN for marks without a track."""
    reg = regions[["chrom", "start", "end"]].drop_duplicates().reset_index(drop=True)
    with tempfile.TemporaryDirectory(prefix="cislink_") as tmp:
        bed = os.path.join(tmp, "regions.bed")   # written once, read by every track
        reg.assign(i=np.arange(len(reg)))[["chrom", "start", "end", "i"]].to_csv(bed, sep="\t", header=False, index=False)

        def run(m):
            p = tracks[m]
            return m, (_sums_bigwig(reg, p) if _is_bigwig(p) else _sums_bedgraph(len(reg), bed, p, bedtools, tmp, m))
        with ThreadPoolExecutor(max_workers=threads) as pool:
            res = dict(pool.map(run, [m for m in MARKS if m in tracks]))
    out = reg.copy()
    for m in MARKS:
        out[m] = np.round(np.asarray(res[m], float), 6) if m in res else np.nan
    return out


def _lookup(sums: pd.DataFrame, chrom, start, end, mark):
    key = pd.MultiIndex.from_frame(sums[["chrom", "start", "end"]])
    return pd.Series(sums[mark].to_numpy(), index=key).reindex(pd.MultiIndex.from_arrays([chrom, start, end])).to_numpy()


# ---------------------------------------------------------------- expression
def _read_tpm(path, id_column):
    t = pd.read_csv(path, sep="\t", usecols=[id_column, "TPM"])
    return t[~t[id_column].astype(str).str.contains("_PAR_Y")]


def transcript_expression(tx_ids, path) -> np.ndarray:
    t = _read_tpm(path, "transcript_id")
    tpm = t.assign(k=t.transcript_id.str.split(".").str[0]).drop_duplicates("k", keep="last").set_index("k").TPM
    return pd.Series(tx_ids).str.split(".").str[0].map(tpm).fillna(0.0).to_numpy()


def gene_expression(genes_symbols, genes: pd.DataFrame, path) -> np.ndarray:
    """Gene TPM via symbol -> versioned gene ID (the last GTF gene record of a symbol)."""
    sym_to_id = dict(zip(genes.gene, genes.gene_id))
    tpm = _read_tpm(path, "gene_id").drop_duplicates("gene_id", keep="last").set_index("gene_id").TPM
    return pd.Series(genes_symbols).map(sym_to_id).map(tpm).fillna(0.0).to_numpy()


def expression_factors(gtf_transcripts: pd.DataFrame, genes: pd.DataFrame, transcript_quant, gene_quant) -> dict:
    """Scale factors that make expression protein-coding TPM (1e6 / protein-coding TPM sum)."""
    pc_tx = set(gtf_transcripts[gtf_transcripts.transcript_type == "protein_coding"].TX.str.split(".").str[0])
    pc_g = set(genes[genes.gene_type == "protein_coding"].gene_id.str.split(".").str[0])
    g = pd.read_csv(gene_quant, sep="\t", usecols=["gene_id", "TPM"]); t = pd.read_csv(transcript_quant, sep="\t", usecols=["transcript_id", "TPM"])
    return hz.expression_factors(g[g.gene_id.str.split(".").str[0].isin(pc_g)].TPM.sum(),
                                 t[t.transcript_id.str.split(".").str[0].isin(pc_tx)].TPM.sum())


# ---------------------------------------------------------------- Hi-C
def _hic_names(hic):
    n = {c.name: c.length for c in hic.getChromosomes()}
    plain = "chr1" not in n and "1" in n
    return (lambda c: c.replace("chr", "") if plain else c), n


def hic_values(hic_path, chrom, a_pos, b_pos, resolution: int = HIC_BIN, norm: str = "KR", matrix_type: str = "oe") -> np.ndarray:
    """Contact between the bin of ``a_pos`` and the bin of ``b_pos`` (same chromosome); NaN where the map has none."""
    import hicstraw
    chrom, a, b = np.asarray(chrom), np.asarray(a_pos, int) // resolution * resolution, np.asarray(b_pos, int) // resolution * resolution
    hic = hicstraw.HiCFile(str(hic_path)); name, lengths = _hic_names(hic)
    want: dict = {}
    for c, e, t in zip(chrom, a, b):
        want.setdefault(c, {}).setdefault(int(e), set()).add(int(t))
    found = {}
    for c, by_e in want.items():
        if name(c) not in lengths:
            continue
        mzd = hic.getMatrixZoomData(name(c), name(c), matrix_type, norm, "BP", resolution)
        for e, ts in by_e.items():
            lo, hi = min(ts | {e}), max(ts | {e})
            for x1, x2, y1, y2 in ((e, e, e, hi), (lo, e, e, e)):
                if x2 < x1 or y2 < y1:
                    continue
                for r in mzd.getRecords(x1, x2, y1, y2):
                    x, y = sorted((r.binX, r.binY))
                    if (x == e and y in ts) or (y == e and x in ts):
                        found[(c, x, y)] = r.counts
    return np.array([found.get((c, min(e, t), max(e, t)), np.nan) for c, e, t in zip(chrom, a, b)], float)


def hic_cell_quantiles(hic_path, resolution: int = HIC_BIN, norm: str = "KR", windows: int = 20, seed: int = 0) -> pd.DataFrame:
    """A map's O/E quantiles per distance stratum, from ``windows`` random 2-Mb windows per autosome."""
    import hicstraw
    hic = hicstraw.HiCFile(str(hic_path)); name, lengths = _hic_names(hic); rng = np.random.default_rng(seed)
    vals, dists = [], []
    for c in AUTOSOMES:
        L, w = lengths[name(c)], 2_000_000
        mzd = hic.getMatrixZoomData(name(c), name(c), "oe", norm, "BP", resolution)
        for s0 in rng.integers(0, max(L - w, 1), windows):
            s0 = int(s0) // resolution * resolution
            m = np.nan_to_num(mzd.getRecordsAsMatrix(s0, s0 + w - resolution, s0, s0 + w - resolution))
            if m.size == 0:
                continue
            idx = np.where(m.sum(0) > 0)[0]
            i, j = np.triu_indices(len(idx)); d = (idx[j] - idx[i]) * resolution; ok = d <= 1e6
            vals.append(m[np.ix_(idx, idx)][i[ok], j[ok]]); dists.append(d[ok])
    return _sig6(hz.hic_quantiles(np.concatenate(vals), np.concatenate(dists)))


# ---------------------------------------------------------------- candidate elements (for chromatin quantiles)
def make_candidates(accessibility_pvalue_bedgraph, out, bedtools: str = "bedtools") -> str:
    """Candidate elements from an accessibility -log10 P bedGraph: bins with -log10 P >= T merged across gaps of
    <= 100 bp, with T the threshold giving closest to 150,000 regions; activity = maximum -log10 P. Writes a BED
    (chrom, start, end, activity) and returns its path."""
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out))) as tmp:
        low = os.path.join(tmp, "low.bedGraph")
        subprocess.run(f"awk -v OFS='\\t' '$4>={min(CANDIDATE_THRESHOLDS)} && $1!~/_/' '{accessibility_pvalue_bedgraph}' > '{low}'", shell=True, check=True)
        best = None
        for t in CANDIDATE_THRESHOLDS:
            p = os.path.join(tmp, f"t{t}.bed")
            subprocess.run(f"awk '$4>={t}' '{low}' | '{bedtools}' merge -i - -d {CANDIDATE_MERGE_GAP} -c 4 -o max > '{p}'", shell=True, check=True)
            n = sum(1 for _ in open(p))
            if best is None or abs(n - CANDIDATE_TARGET) < abs(best[1] - CANDIDATE_TARGET):
                best = (t, n, p)
        os.replace(best[2], out)
    return str(out)


def cell_chromatin_quantiles(candidates_bed, tracks: dict[str, str], bedtools: str = "bedtools", threads: int = 4) -> pd.DataFrame:
    """Quantiles of each track's per-bp mean signal over the cell's candidate elements."""
    cand = pd.read_csv(candidates_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    reg = cand.drop_duplicates().reset_index(drop=True)
    sums = signal_sums(reg, tracks, bedtools, threads)
    width = (reg.end - reg.start).clip(lower=1).to_numpy()
    return _sig6(hz.chromatin_quantiles({m: sums[m].to_numpy() / width for m in MARKS if m in tracks}))


# ---------------------------------------------------------------- distance and gene context
def distance_features(rows: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """distance, log10_distance, Enh_length, powerlaw_contact, n_genes_between, gene_rank, is_nearest_gene and
    ubiquitous_gene for rows with chrom, es, ee, ps, pe, gene."""
    mid, tss = (rows.es + rows.ee) / 2, (rows.ps + rows.pe) / 2
    dist = (mid - tss).abs().clip(lower=1)
    out = pd.DataFrame({"distance": dist, "log10_distance": np.log10(dist), "Enh_length": (rows.ee - rows.es).clip(lower=1),
                        "powerlaw_contact": (np.maximum(dist, HIC_BIN) / HIC_BIN) ** -POWERLAW_GAMMA}, index=rows.index)
    g = genes[genes.gene_type == "protein_coding"].sort_values(["chrom", "tss"])
    by = {c: (d.tss.to_numpy(), d.gene.to_numpy()) for c, d in g.groupby("chrom")}
    between, rank = np.zeros(len(rows), int), np.zeros(len(rows), int)
    for i, (c, m, t, gn) in enumerate(zip(rows.chrom, mid, tss, rows.gene)):
        pos, names = by.get(c, (np.array([]), np.array([])))
        other = names != gn
        lo, hi = min(m, t), max(m, t)
        between[i] = int(((pos > lo) & (pos < hi) & other).sum())
        rank[i] = 1 + int((np.abs(pos - m)[other] < abs(m - t)).sum())
    out["n_genes_between"], out["gene_rank"], out["is_nearest_gene"] = between, rank, (rank == 1).astype(int)
    out["ubiquitous_gene"] = rows.gene.isin(ubiquitous_genes()).astype(int).to_numpy()
    return out


# ---------------------------------------------------------------- full build
def build_features(rows: pd.DataFrame, *, gtf, tracks: dict[str, str], transcript_quantification, gene_quantification,
                   candidates_bed, hic=None, hic_resolution: int = HIC_BIN, hic_norm: str = "KR", map_hic: bool = True,
                   megamap=MEGAMAP, bedtools: str = "bedtools", threads: int = 4, cell_quantiles: pd.DataFrame | None = None,
                   hic_quantiles: pd.DataFrame | None = None, log=print) -> pd.DataFrame:
    """Harmonized feature table for ``rows`` (from ``variant_pairs`` or ``element_pairs``), ready for scoring.

    ``hic``: the cell's .hic for CisLink-E2G (omit for V2G only). ``map_hic``: map the cell's O/E onto K562 within
    distance strata (default; the manuscript did this for the lab's 2-2.5-kb maps). ``megamap``: .hic path or URL of
    the multi-tissue megamap for CisLink-V2G (``None`` to skip)."""
    transcripts, genes = read_gtf(gtf)
    out = rows[["ID", "chrom", "es", "ee", "ps", "pe", "gene", "TX"]].copy()
    out["Enh"] = out.chrom + ":" + out.es.astype(str) + ":" + out.ee.astype(str)
    out["Prom"] = out.chrom + ":" + out.ps.astype(str) + ":" + out.pe.astype(str)
    log(f"{len(out)} rows; signals over {len(tracks)} tracks")
    regions = pd.concat([out[["chrom", "es", "ee"]].set_axis(["chrom", "start", "end"], axis=1),
                         out[["chrom", "ps", "pe"]].set_axis(["chrom", "start", "end"], axis=1)])
    sums = signal_sums(regions, tracks, bedtools, threads)
    for m in MARKS:
        out[f"{m}_Enh"] = _lookup(sums, out.chrom, out.es, out.ee, m)
        out[f"{m}_Prom"] = _lookup(sums, out.chrom, out.ps, out.pe, m)
    out["Tx_expression"] = transcript_expression(out.TX, transcript_quantification)
    out["Gene_expression"] = gene_expression(out.gene, genes, gene_quantification)
    # harmonize chromatin and expression
    if cell_quantiles is None:
        log("cell chromatin quantiles over candidate elements")
        cell_quantiles = cell_chromatin_quantiles(candidates_bed, tracks, bedtools, threads)
    out = hz.harmonize_chromatin(out, cell_quantiles, enh_len=out.ee - out.es, prom_len=out.pe - out.ps)
    out = hz.harmonize_expression(out, expression_factors(transcripts, genes, transcript_quantification, gene_quantification))
    out = pd.concat([out, distance_features(out, genes)], axis=1)
    emid, tss = (out.es + out.ee) // 2, (out.ps + out.pe) // 2
    if hic is not None:
        log(f"cell Hi-C O/E at {hic_resolution} bp")
        v = hic_values(hic, out.chrom, emid, tss, hic_resolution, hic_norm, "oe")
        if map_hic:
            hq = hic_quantiles if hic_quantiles is not None else hic_cell_quantiles(hic, hic_resolution, hic_norm)
            bd = np.abs(emid // hic_resolution * hic_resolution - tss // hic_resolution * hic_resolution)
            v = hz.harmonize_hic(v, bd, hq)
        out["HiC_Contact"] = v
    else:
        out["HiC_Contact"] = np.nan
    if megamap is not None:
        log("megamap O/E at 5 kb")
        out["HiC_megamap_oe"] = hic_values(megamap, out.chrom, emid, tss, HIC_BIN, "SCALE", "oe")
    else:
        out["HiC_megamap_oe"] = np.nan   # scored as no contact
    return out
