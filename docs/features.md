# Input features

Both models take one row per element–gene (or variant–gene) pair with the columns below. Signals are summed over the
element (or the ±250-bp variant window) and over a ±250-bp promoter window around each GENCODE v29 TSS, then
harmonized to the K562 scale (`cislink.harmonize`). E2G uses the 40 features with an E2G coefficient; V2G uses the 16
with a V2G coefficient. Accessibility comes from ATAC-seq or DNase-seq; the earlier column names `DNase_Enh` and
`DNase_Prom` are still accepted.

Three features are new in v1.3:

- `F1_dnase` and `F1_h3k27ac` (cross-epigenome correlation) do not depend on the scored cell. They are computed from the
  Roadmap activity reference (`cislink reference download`) by `cislink build` or `cislink add-f1`
  (`cislink.roadmap`). The epigenomes are those with observed, not imputed, DNase-seq or H3K27ac and RNA-seq, minus
  E000 and E123 (K562, the training cell). Imputed data are not used because imputed signals are predicted from the
  same epigenome's other marks, which would make the correlation partly circular.
- `DxA` is computed by the package from `log10_distance` and `Accessibility_Enh` when a table lacks it.

Missing Hi-C contact and missing F1 (an element that does not lift to hg19, or a gene without Roadmap RNA-seq) are
treated as 0, as in training; any other missing value is an error.

| Column | Label | Group | Definition | Harmonization | E2G coef. | V2G coef. |
|---|---|---|---|---|---|---|
| `Accessibility_Enh` | Accessibility (element) | accessibility | accessibility (DNase-seq or ATAC-seq) summed over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.346 | +0.617 |
| `Accessibility_Prom` | Accessibility (promoter) | accessibility | accessibility over the ±250-bp promoter window | quantile-mapped to K562 over genome-wide candidate elements | -0.129 | -0.459 |
| `Tx_expression` | Transcript expression | expression | expression (TPM over protein-coding transcripts) of the gene's transcript | recomputed over protein-coding transcripts | +0.112 | +0.117 |
| `Gene_expression` | Gene expression | expression | gene expression (TPM, protein-coding) | recomputed over protein-coding transcripts | -0.254 | -0.292 |
| `HiC_megamap_oe` | Hi-C contact (megamap) | contact | observed/expected contact in the multi-tissue megamap, 5 kb | none |  | +0.193 |
| `distance` | Distance | distance and gene context | element–TSS distance (bp) | none | -0.168 | -0.164 |
| `log10_distance` | Log distance | distance and gene context | log10 element–TSS distance | none | -0.554 | -0.383 |
| `n_genes_between` | Genes between element and TSS | distance and gene context | number of genes between element and TSS | none | -0.839 | -0.979 |
| `gene_rank` | Distance rank of the gene | distance and gene context | rank of the gene by distance among genes near the element | none | -0.212 | -0.226 |
| `is_nearest_gene` | Nearest gene | distance and gene context | gene is the nearest to the element | none | +0.127 | +0.118 |
| `ubiquitous_gene` | Ubiquitously expressed gene | expression | gene is ubiquitously expressed | none | -0.545 | -0.582 |
| `Enh_length` | Element length | element | element length | none | -0.428 | -0.463 |
| `powerlaw_contact` | Power-law contact | contact | power-law contact expected from distance | none | +0.044 | +0.130 |
| `F1_dnase` | Accessibility–expression correlation across epigenomes | cross-epigenome correlation | Spearman correlation, across 28 Roadmap epigenomes (K562 excluded), between the element's DNase-seq signal and the gene's expression | none (fixed Roadmap reference, independent of the scored cell) | +0.108 | +0.149 |
| `F1_h3k27ac` | H3K27ac–expression correlation across epigenomes | cross-epigenome correlation | Spearman correlation, across 47 Roadmap epigenomes (K562 excluded), between the element's H3K27ac signal and the gene's expression | none (fixed Roadmap reference, independent of the scored cell) | +0.161 | +0.107 |
| `DxA` | Log distance × accessibility | accessibility | log10 element–TSS distance × element accessibility | product of the quantile-mapped accessibility and log10 distance | +0.286 | +0.579 |
| `H4K20me1_Enh` | H4K20me1 (element) | histone marks | H4K20me1 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.386 |  |
| `H4K20me1_Prom` | H4K20me1 (promoter) | histone marks | H4K20me1 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.097 |  |
| `H3K9me3_Enh` | H3K9me3 (element) | histone marks | H3K9me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.029 |  |
| `H3K9me3_Prom` | H3K9me3 (promoter) | histone marks | H3K9me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.032 |  |
| `H3K9ac_Enh` | H3K9ac (element) | histone marks | H3K9ac ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.414 |  |
| `H3K9ac_Prom` | H3K9ac (promoter) | histone marks | H3K9ac ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.065 |  |
| `H3K79me2_Enh` | H3K79me2 (element) | histone marks | H3K79me2 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.088 |  |
| `H3K79me2_Prom` | H3K79me2 (promoter) | histone marks | H3K79me2 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.022 |  |
| `H3K4me3_Enh` | H3K4me3 (element) | histone marks | H3K4me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.637 |  |
| `H3K4me3_Prom` | H3K4me3 (promoter) | histone marks | H3K4me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.157 |  |
| `H3K4me2_Enh` | H3K4me2 (element) | histone marks | H3K4me2 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.508 |  |
| `H3K4me2_Prom` | H3K4me2 (promoter) | histone marks | H3K4me2 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.004 |  |
| `H3K4me1_Enh` | H3K4me1 (element) | histone marks | H3K4me1 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.308 |  |
| `H3K4me1_Prom` | H3K4me1 (promoter) | histone marks | H3K4me1 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.068 |  |
| `H3K36me3_Enh` | H3K36me3 (element) | histone marks | H3K36me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.069 |  |
| `H3K36me3_Prom` | H3K36me3 (promoter) | histone marks | H3K36me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.025 |  |
| `H3K27me3_Enh` | H3K27me3 (element) | histone marks | H3K27me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.086 |  |
| `H3K27me3_Prom` | H3K27me3 (promoter) | histone marks | H3K27me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.096 |  |
| `H3K27ac_Enh` | H3K27ac (element) | histone marks | H3K27ac ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +1.428 |  |
| `H3K27ac_Prom` | H3K27ac (promoter) | histone marks | H3K27ac ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.467 |  |
| `H2A.Z_Enh` | H2A.Z (element) | histone marks | H2A.Z ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.530 |  |
| `H2A.Z_Prom` | H2A.Z (promoter) | histone marks | H2A.Z ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.029 |  |
| `CTCF_Enh` | CTCF (element) | CTCF | CTCF ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.136 |  |
| `CTCF_Prom` | CTCF (promoter) | CTCF | CTCF ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.179 |  |
| `HiC_Contact` | Hi-C contact (cell's own map) | contact | observed/expected contact in the cell's own Hi-C, mapped to K562 5 kb within distance strata | distance-stratified quantile mapping to K562 5 kb | +0.472 |  |

Coefficients are standardized (after the `sign(x)·log1p(|x|)` transform and standardization).
