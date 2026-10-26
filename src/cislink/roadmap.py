"""Cross-epigenome activity-expression correlation features of CisLink v1.3 (``F1_dnase``, ``F1_h3k27ac``).

For an element-gene pair, F1 is the Spearman correlation, across Roadmap Epigenomics reference epigenomes, between the
element's signal and the gene's expression:

* element signal: the integral of the Roadmap consolidated MACS2 -log10 P signal over the element (hg19), DNase-seq for
  ``F1_dnase`` (28 epigenomes) and H3K27ac for ``F1_h3k27ac`` (47 epigenomes);
* gene expression: log(1 + RPKM) of the Roadmap 57-epigenome protein-coding RNA-seq matrix;
* epigenomes: those with observed (not imputed) data for the mark and RNA-seq, excluding E000 (a pooled reference RNA)
  and E123 (K562, the training cell).

The values come from the Roadmap activity reference (``cislink reference download``; 16.3 GB): every track binned at
25 bp on hg19; an element's signal is the sum of its bins, edge bins prorated by the covered fraction. GRCh38 elements are
lifted to hg19 with CrossMap and the UCSC hg38ToHg19 chain (one-to-one only). A pair gets 0 when its element does not lift,
its gene has no RNA-seq entry, or either vector is constant. F1 does not depend on the user's cell: cell-specific evidence
comes from the cell's own inputs.
"""
from __future__ import annotations

import gzip
import os
import re
import shlex
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

ENV = "CISLINK_ROADMAP_REF"
MARKS = (("DNase", "F1_dnase"), ("H3K27ac", "F1_h3k27ac"))
F1_FEATURES = tuple(f for _, f in MARKS)
CHAIN = "hg38ToHg19.over.chain.gz"


class ReferenceError(FileNotFoundError):
    pass


def reference_path(path=None) -> str:
    """The reference directory: ``path``, else ``$CISLINK_ROADMAP_REF``; a clear error when neither holds it."""
    p = path or os.environ.get(ENV)
    if not p:
        raise ReferenceError(f"the CisLink v1.3 models need the Roadmap activity reference: run 'cislink reference download "
                             f"--dest DIR' and pass --roadmap-ref DIR (or set {ENV}=DIR, or 'roadmap_reference: DIR' in the config)")
    if not os.path.isfile(os.path.join(p, "manifest.tsv")):
        raise ReferenceError(f"{p} is not a Roadmap activity reference (no manifest.tsv); see 'cislink reference download'")
    return str(p)


def _parse(el: str):
    c, s, e = str(el).split(":")
    return c, int(s), int(e)


def _z(r: pd.DataFrame) -> pd.DataFrame:
    r = r.sub(r.mean(axis=1), axis=0)
    return r.div(np.sqrt((r ** 2).sum(axis=1)), axis=0)


class RoadmapActivity:
    """Reader of the Roadmap activity reference."""

    def __init__(self, path=None):
        self.ref = reference_path(path)
        self.man = pd.read_csv(os.path.join(self.ref, "manifest.tsv"), sep="\t")
        self.size = dict(zip(self.man.track, self.man.bin_bp.astype(int)))
        self.rpkm = pd.read_csv(os.path.join(self.ref, "rna_rpkm_55epigenomes.tsv.gz"), sep="\t", index_col=0)
        self.eids = {m: list(self.man.loc[self.man.mark == m, "epigenome"]) for m, _ in MARKS}

    def _track_values(self, track: str, elements) -> np.ndarray:
        """Integral of -log10 P over each (chrom, start, end), edge bins prorated by their covered fraction."""
        z = np.load(os.path.join(self.ref, "tracks", f"{track}.npz")); out = np.zeros(len(elements))
        by: dict = {}
        for i, (c, _, _) in enumerate(elements):
            by.setdefault(c, []).append(i)
        size = self.size[track]
        for c, idx in by.items():
            if c not in z.files:
                continue
            b = z[c]; cs = np.concatenate([[0.0], np.cumsum(b, dtype=np.float64)]); n = len(b)
            for i in idx:
                _, s, e = elements[i]
                if e <= s or s < 0:
                    continue
                b0, b1 = s // size, (e - 1) // size
                if b0 >= n:
                    continue
                b1 = min(b1, n - 1)
                if b0 == b1:
                    out[i] = b[b0] * (e - s) / size
                else:
                    out[i] = cs[b1] - cs[b0 + 1] + b[b0] * ((b0 + 1) * size - s) / size + b[b1] * (e - b1 * size) / size
        return out

    def signals(self, elements_hg19, mark: str) -> pd.DataFrame:
        """Elements x epigenomes signal for one mark ('DNase' or 'H3K27ac'); elements as hg19 'chr:start:end'."""
        els = [_parse(x) for x in elements_hg19]
        return pd.DataFrame({e: self._track_values(f"{e}-{mark}", els) for e in self.eids[mark]}, index=list(elements_hg19))

    def f1_table(self, elements_hg19, gene_ids) -> pd.DataFrame:
        """F1_dnase and F1_h3k27ac for aligned hg19 elements ('chr:start:end', NaN if unmapped) and stable Ensembl gene IDs
        (NaN if unknown); missing -> 0."""
        el = pd.Series(list(elements_hg19), dtype=object); gid = pd.Series(list(gene_ids), dtype=object)
        uniq = sorted(set(el.dropna()))
        out = pd.DataFrame(index=range(len(el)))
        for mark, name in MARKS:
            eids = self.eids[mark]
            S = _z(self.signals(uniq, mark).rank(axis=1)) if uniq else pd.DataFrame(columns=eids)
            X = _z(np.log1p(self.rpkm[eids]).rank(axis=1))
            ok = el.isin(S.index) & gid.isin(X.index); v = np.zeros(len(el))
            if ok.any():
                v[ok.values] = np.nansum(S.loc[el[ok]].values * X.loc[gid[ok]].values, axis=1)
            v[~np.isfinite(v)] = 0.0
            out[name] = v
        return out


# ---------------------------------------------------------------- GRCh38 -> hg19 and gene IDs
def lift_to_hg19(elements, chain: str, crossmap: str = "CrossMap", tmpdir=None) -> dict:
    """GRCh38 'chr:start:end' -> hg19 'chr:start:end' with CrossMap (bed); unmapped or split elements -> NaN."""
    u = pd.Series(pd.unique(pd.Series(list(elements), dtype=object).dropna()))
    if not len(u):
        return {}
    exe = shutil.which(crossmap)
    if not exe:
        raise RuntimeError(f"{crossmap} not found; it lifts elements to hg19 for the F1 features "
                           "(pixi install -e macs2 provides .pixi/envs/macs2/bin/CrossMap; or set 'crossmap:' in the config)")
    p = u.str.split(":", expand=True)
    with tempfile.TemporaryDirectory(prefix="cislink_lift_", dir=tmpdir) as t:
        a, b = os.path.join(t, "a.bed"), os.path.join(t, "b.bed")
        pd.DataFrame({"c": p[0], "s": p[1].astype(int), "e": p[2].astype(int), "n": np.arange(len(u))}).to_csv(a, sep="\t", header=False, index=False)
        subprocess.run(f"{shlex.quote(exe)} bed {shlex.quote(chain)} {shlex.quote(a)} {shlex.quote(b)} > {shlex.quote(os.path.join(t, 'log'))} 2>&1",
                       shell=True, check=True)
        if os.path.getsize(b):
            r = pd.read_csv(b, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["c", "s", "e", "n"])
        else:
            r = pd.DataFrame(columns=["c", "s", "e", "n"])
    r = r[~r.n.duplicated(keep=False)]   # one-to-one only
    m = dict(zip(r.n, r.c.astype(str) + ":" + r.s.astype(str) + ":" + r.e.astype(str)))
    return dict(zip(u, [m.get(i, np.nan) for i in range(len(u))]))


def gene_ids(gtf, transcripts, symbols) -> np.ndarray:
    """Stable Ensembl gene ID of each row: from the transcript (GTF transcript -> gene), else from the symbol (a
    protein-coding gene preferred)."""
    tx2g, sym = {}, {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            g = re.search(r'gene_id "([^".]+)', f[8])
            if f[2] == "transcript":
                tx2g[re.search(r'transcript_id "([^".]+)', f[8]).group(1)] = g.group(1)
            elif f[2] == "gene":
                name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
                pc = 'gene_type "protein_coding"' in f[8]
                if name not in sym or (pc and not sym[name][1]):
                    sym[name] = (g.group(1), pc)
    t = pd.Series(list(transcripts), dtype=object).astype(str).str.split(".").str[0].map(tx2g)
    s = pd.Series(list(symbols), dtype=object).map(lambda x: sym.get(x, (np.nan,))[0])
    return t.fillna(s).to_numpy()


def chain_path(ref: str) -> str:
    p = os.path.join(ref, CHAIN)
    if not os.path.isfile(p):
        raise ReferenceError(f"{p} missing; the reference ships the UCSC hg38ToHg19 chain (cislink reference download)")
    return p


def add_f1(table: pd.DataFrame, gtf, ref=None, crossmap: str = "CrossMap", genome: str = "GRCh38", log=print) -> pd.DataFrame:
    """Copy of a feature table (columns Enh 'chr:start:end', TX and, optionally, gene for the symbol fallback) with
    F1_dnase and F1_h3k27ac added."""
    ra = RoadmapActivity(ref)
    out = table.copy()
    if genome == "hg19":
        el = out.Enh
    else:
        log("F1: lifting elements to hg19 (CrossMap)")
        el = out.Enh.map(lift_to_hg19(out.Enh, chain_path(ra.ref), crossmap))
    log(f"F1: {el.notna().mean():.1%} of rows lifted; correlations across {len(ra.eids['DNase'])} DNase and "
        f"{len(ra.eids['H3K27ac'])} H3K27ac epigenomes")
    f = ra.f1_table(el.to_numpy(), gene_ids(gtf, out.TX, out.gene if "gene" in out else [np.nan] * len(out)))
    out["F1_dnase"], out["F1_h3k27ac"] = f.F1_dnase.to_numpy(), f.F1_h3k27ac.to_numpy()
    return out


# ---------------------------------------------------------------- download
def release_info() -> dict:
    """Archive names, sizes and md5 of the published reference (src/cislink/reference/roadmap_reference.json)."""
    import json
    from importlib import resources
    with resources.files("cislink.reference").joinpath("roadmap_reference.json").open() as fh:
        return json.load(fh)


def _md5(path, block: int = 1 << 24) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def verify(ref) -> int:
    """Check every file of an extracted reference against files_md5.tsv (bytes and md5); returns the number checked."""
    t = pd.read_csv(os.path.join(ref, "files_md5.tsv"), sep="\t")
    bad = []
    for f, n, h in zip(t.file, t.bytes, t.md5):
        p = os.path.join(ref, f)
        if not os.path.isfile(p) or os.path.getsize(p) != int(n) or _md5(p) != h:
            bad.append(f)
    if bad:
        raise ReferenceError(f"{len(bad)} reference files missing or corrupt, e.g. {bad[:3]}; download again")
    return len(t)


def download(dest, base_url: str | None = None, keep_archives: bool = False, log=print) -> str:
    """Download the Roadmap activity reference archives into ``dest``, check their md5, extract and verify every file.
    Needs about 16.3 GB (plus the archives while extracting). Returns the reference directory."""
    import tarfile
    import urllib.request
    info = release_info()
    base = base_url or info.get("base_url")
    if not base:
        raise ReferenceError(f"the reference is not published yet (DOI {info.get('doi')}); pass --base-url with the location of the archives")
    os.makedirs(dest, exist_ok=True)
    for a in info["archives"]:
        path = os.path.join(dest, a["name"])
        if not (os.path.isfile(path) and os.path.getsize(path) == a["bytes"] and _md5(path) == a["md5"]):
            src = base.rstrip("/") + "/" + a["name"]
            log(f"downloading {a['name']} ({a['bytes'] / 1e9:.1f} GB)")
            if src.startswith(("http://", "https://")):
                urllib.request.urlretrieve(src + info.get("url_suffix", ""), path)
            else:
                shutil.copyfile(src, path)
            if _md5(path) != a["md5"]:
                raise ReferenceError(f"{a['name']}: md5 mismatch after download")
        log(f"extracting {a['name']}")
        with tarfile.open(path) as tf:
            if hasattr(tarfile, "data_filter"):
                tf.extractall(dest, filter="data")
            else:   # Python < 3.10.12: refuse members outside dest
                root = os.path.realpath(dest)
                for mem in tf.getmembers():
                    if not os.path.realpath(os.path.join(dest, mem.name)).startswith(root + os.sep) or mem.issym() or mem.islnk():
                        raise ReferenceError(f"unsafe archive member {mem.name}")
                tf.extractall(dest)
        if not keep_archives:
            os.remove(path)
    ref = os.path.join(dest, info["directory"])
    log(f"verified {verify(ref)} files in {ref}")
    return ref
