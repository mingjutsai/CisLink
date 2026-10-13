# CisLink

[![CI](https://github.com/mingjutsai/CisLink/actions/workflows/ci.yml/badge.svg)](https://github.com/mingjutsai/CisLink/actions/workflows/ci.yml)

**General models that link noncoding variants and cis-regulatory elements to their target genes.**

CisLink has two models, both trained once on 10,250 CRISPR-tested element–gene pairs in K562 cells and applied
across cell types, traits and accessibility assays:

| | **CisLink-E2G** | **CisLink-V2G** |
|---|---|---|
| Question | Which elements regulate which genes in this cell? | Which gene does this GWAS variant regulate? |
| Inputs (features) | accessibility, 11 histone marks, CTCF, expression, the cell's own Hi-C, distance and gene context (37) | accessibility, expression, multi-tissue Hi-C megamap, distance and gene context (13) |
| Needs | ATAC-seq or DNase-seq, histone ChIP-seq, RNA-seq, Hi-C | ATAC-seq or DNase-seq, RNA-seq |
| Use it to | map enhancer–gene links, especially for elements 10–200 kb from the gene; rank elements for CRISPR follow-up | score every fine-mapped variant at a GWAS locus, in a peak or not, and call a gene at every locus |

Both are L2-regularized logistic regressions stored as plain JSON (`src/cislink/models/`), so scoring needs only NumPy.

## Install

With [pixi](https://pixi.sh) (recommended; reproduces the tested environment):

```bash
git clone https://github.com/mingjutsai/CisLink.git
cd CisLink
pixi install          # Python, NumPy, pandas, scikit-learn, Jupyter
pixi run test         # unit tests
```

With pip:

```bash
pip install git+https://github.com/mingjutsai/CisLink.git
```

## Tutorial

### 1. Score a feature table

A feature table has one row per element–gene (or variant–gene) pair and the columns listed in
[`docs/features.md`](docs/features.md); `examples/example_features.tsv` is a small example.

```python
import pandas as pd
from cislink import load_model, score

x = pd.read_csv("examples/example_features.tsv", sep="\t")
s = score(x)                                  # adds CisLink-E2G and CisLink-V2G columns
e2g = load_model("e2g")
s["E2G link"] = s["CisLink-E2G"] >= e2g.threshold(70)
```

From the shell:

```bash
pixi run cislink info                                                   # features and thresholds
pixi run cislink score examples/example_features.tsv --keep ID --threshold 70 -o scores.tsv
```

When a gene has several transcripts (several rows per pair), take the maximum over its rows.

### 2. Read the scores

Scores are probabilities on the scale of the K562 CRISPR training pairs (4.6% positive); they are not absolute
probabilities in a new cell type. Use them to **rank** candidates, or against the **recall-based thresholds**, the
score at which held-out K562 CRISPR pairs reach a given recall:

| | 70% recall | 50% recall |
|---|---|---|
| CisLink-E2G | 0.180 | 0.529 |
| CisLink-V2G | 0.139 | 0.410 |

The scale transfers to HCT116 and GM12878 but scores are overconfident in osteoblasts (notebook 04); the 70%-recall
threshold kept 64–97% recall across cell types. `model.contributions(table)` gives each feature's contribution to the
log-odds of a pair.

### 3. Build features for a new cell type

`cislink build` computes the harmonized feature table from a cell type's raw data (GTF, bedGraph or bigWig tracks,
RNA-seq quantifications, `.hic`), scores it, and `cislink rank` ranks genes per variant or element. It reproduces the
manuscript's matrices exactly; see [`docs/pipeline.md`](docs/pipeline.md) and `examples/config.example.yaml`.

```bash
pixi install -e pipeline
pixi run -e pipeline cislink validate --variants variants.tsv --config my_cell.yaml   # preflight
pixi run -e pipeline cislink build --variants variants.tsv --config my_cell.yaml --score -o features.tsv
pixi run cislink rank features.tsv -o genes.tsv
```

The harmonization steps are also available on their own (`cislink.harmonize`), with the K562 references shipped
in the package:

```python
from cislink import harmonize as hz

# chromatin: per-bp mean signal over the cell's genome-wide candidate elements -> 1,001 quantiles per track
cell_q = hz.chromatin_quantiles({"DNase": dnase_means, "H3K27ac": h3k27ac_means})
x = hz.harmonize_chromatin(x, cell_q, enh_len=enh_lengths, prom_len=prom_lengths)

# expression: protein-coding TPM
x = hz.harmonize_expression(x, hz.expression_factors(gene_pc_tpm_sum, transcript_pc_tpm_sum))

# Hi-C (CisLink-E2G): observed/expected, mapped onto K562 5 kb within 16 distance strata
hic_q = hz.hic_quantiles(sampled_oe, sampled_bin_distance)
x["HiC_Contact"] = hz.harmonize_hic(x.HiC_Contact, bin_distance, hic_q)
```

Distance and gene-context features are defined in `docs/features.md`. CisLink-V2G uses the multi-tissue Hi-C
megamap recommended by ENCODE-rE2G for cell types without Hi-C (`HiC_megamap_oe`).

## Reproducible analyses

The notebooks run on the data in `data/` (see [`data/README.md`](data/README.md)):

| Notebook | Reproduces |
|---|---|
| [`01_quickstart`](notebooks/01_quickstart.ipynb) | scoring, thresholds, per-feature contributions, CLI |
| [`02_crispr_benchmark`](notebooks/02_crispr_benchmark.ipynb) | Table 1a (K562, HCT116, GM12878), paired differences with gene-clustered CIs, Table 1c (by element–TSS distance), pooled comparison |
| [`03_feature_importance`](notebooks/03_feature_importance.ipynb) | grouped permutation importance, overall and by distance (Fig. 5b,c) |
| [`04_calibration_thresholds`](notebooks/04_calibration_thresholds.ipynb) | calibration across cell types and the recall-based thresholds |
| [`05_gwas_silver_standard`](notebooks/05_gwas_silver_standard.ipynb) | 568 credible sets for 31 traits: top-1 and MRR (Fig. 3a), by trait (3b), locus-level precision–recall without L2G's training genes (3c), nearest versus non-nearest sets (3d), gene–locus pairs against the nearest gene |

```bash
pixi run lab          # open the notebooks
pixi run notebooks    # re-execute all of them
```

Key results: on held-out K562 CRISPR pairs, CisLink-E2G AUPRC 0.688 (ENCODE-rE2G Extended 0.678, CisLink-V2G 0.578,
distance 0.366), and 0.683 against 0.158 for distance 50–200 kb from the TSS; on the 568 silver-standard credible
sets, CisLink-V2G ranks the silver gene first in 52.6% (nearest gene 51.8%, L2G 45.1%, PoPS 37.1%, ENCODE-rE2G 16.5%).

## Development

```bash
pixi run test                     # pytest
pixi run -e pipeline pytest       # including the feature pipeline (bedtools, pyBigWig, hic-straw)
pixi run lint                     # ruff
pixi run -e build build           # wheel and sdist in dist/
```

Continuous integration (GitHub Actions) runs the tests on Linux, macOS and Windows with Python 3.10 and 3.12,
the feature-pipeline tests on Linux and macOS, lints the code, executes every notebook, and builds and installs the wheel. Pushing a tag `v*` runs the tests,
checks the tag against the package version, and attaches the wheel, sdist and model files to a GitHub release.

## Citation

Tsai, M.-J. *et al.* General models link noncoding variants and cis-regulatory elements to their target genes.
(2026, in preparation). See `CITATION.cff`.

## License

Code and models: MIT (see `LICENSE`). Data in `data/` come from the ENCODE CRISPR benchmark and from the silver
standard of Weeks et al. and Gschwind et al.; cite the original studies listed in `data/README.md`.
