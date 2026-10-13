# Input features

Both models take one row per element–gene (or variant–gene) pair with the columns below. Signals are summed over the
element (or the ±250-bp variant window) and over a ±250-bp promoter window around each GENCODE v29 TSS, then
harmonized to the K562 scale (`cislink.harmonize`). E2G uses all 37 features; V2G uses the 13 with a V2G coefficient.
Missing Hi-C contact is treated as 0; any other missing value is an error. Accessibility comes from ATAC-seq or DNase-seq;
the earlier column names `DNase_Enh` and `DNase_Prom` are still accepted.

| Column | Label | Group | Definition | Harmonization | E2G coef. | V2G coef. |
|---|---|---|---|---|---|---|
| `Accessibility_Enh` | Accessibility (element) | accessibility | accessibility (DNase-seq or ATAC-seq) summed over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.620 | +1.212 |
| `Accessibility_Prom` | Accessibility (promoter) | accessibility | accessibility over the ±250-bp promoter window | quantile-mapped to K562 over genome-wide candidate elements | -0.140 | -0.478 |
| `Tx_expression` | Transcript expression | expression | expression (TPM over protein-coding transcripts) of the gene's transcript | recomputed over protein-coding transcripts | +0.153 | +0.156 |
| `Gene_expression` | Gene expression | expression | gene expression (TPM, protein-coding) | recomputed over protein-coding transcripts | -0.266 | -0.304 |
| `HiC_megamap_oe` | Hi-C contact (megamap) | contact | observed/expected contact in the multi-tissue megamap, 5 kb | nan |  | +0.202 |
| `distance` | Distance | distance and gene context | element–TSS distance (bp) | nan | -0.186 | -0.062 |
| `log10_distance` | Log distance | distance and gene context | log10 element–TSS distance | nan | -0.500 | -0.430 |
| `n_genes_between` | Genes between element and TSS | distance and gene context | number of genes between element and TSS | nan | -0.850 | -1.074 |
| `gene_rank` | Distance rank of the gene | distance and gene context | rank of the gene by distance among genes near the element | nan | -0.254 | -0.186 |
| `is_nearest_gene` | Nearest gene | distance and gene context | gene is the nearest to the element | nan | +0.139 | +0.136 |
| `ubiquitous_gene` | Ubiquitously expressed gene | expression | gene is ubiquitously expressed | nan | -0.574 | -0.624 |
| `Enh_length` | Element length | element | element length | nan | -0.414 | -0.466 |
| `powerlaw_contact` | Power-law contact | contact | power-law contact expected from distance | nan | +0.052 | +0.128 |
| `H4K20me1_Enh` | H4K20me1 (element) | histone marks | H4K20me1 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.382 |  |
| `H4K20me1_Prom` | H4K20me1 (promoter) | histone marks | H4K20me1 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.108 |  |
| `H3K9me3_Enh` | H3K9me3 (element) | histone marks | H3K9me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.030 |  |
| `H3K9me3_Prom` | H3K9me3 (promoter) | histone marks | H3K9me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.030 |  |
| `H3K9ac_Enh` | H3K9ac (element) | histone marks | H3K9ac ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.358 |  |
| `H3K9ac_Prom` | H3K9ac (promoter) | histone marks | H3K9ac ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.039 |  |
| `H3K79me2_Enh` | H3K79me2 (element) | histone marks | H3K79me2 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.103 |  |
| `H3K79me2_Prom` | H3K79me2 (promoter) | histone marks | H3K79me2 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.028 |  |
| `H3K4me3_Enh` | H3K4me3 (element) | histone marks | H3K4me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.612 |  |
| `H3K4me3_Prom` | H3K4me3 (promoter) | histone marks | H3K4me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.151 |  |
| `H3K4me2_Enh` | H3K4me2 (element) | histone marks | H3K4me2 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.472 |  |
| `H3K4me2_Prom` | H3K4me2 (promoter) | histone marks | H3K4me2 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.005 |  |
| `H3K4me1_Enh` | H3K4me1 (element) | histone marks | H3K4me1 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.324 |  |
| `H3K4me1_Prom` | H3K4me1 (promoter) | histone marks | H3K4me1 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.084 |  |
| `H3K36me3_Enh` | H3K36me3 (element) | histone marks | H3K36me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.079 |  |
| `H3K36me3_Prom` | H3K36me3 (promoter) | histone marks | H3K36me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.018 |  |
| `H3K27me3_Enh` | H3K27me3 (element) | histone marks | H3K27me3 ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +0.073 |  |
| `H3K27me3_Prom` | H3K27me3 (promoter) | histone marks | H3K27me3 ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.088 |  |
| `H3K27ac_Enh` | H3K27ac (element) | histone marks | H3K27ac ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | +1.358 |  |
| `H3K27ac_Prom` | H3K27ac (promoter) | histone marks | H3K27ac ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.474 |  |
| `H2A.Z_Enh` | H2A.Z (element) | histone marks | H2A.Z ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.531 |  |
| `H2A.Z_Prom` | H2A.Z (promoter) | histone marks | H2A.Z ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | +0.015 |  |
| `CTCF_Enh` | CTCF (element) | CTCF | CTCF ChIP-seq over the element | quantile-mapped to K562 over genome-wide candidate elements | -0.121 |  |
| `CTCF_Prom` | CTCF (promoter) | CTCF | CTCF ChIP-seq over the promoter | quantile-mapped to K562 over genome-wide candidate elements | -0.185 |  |
| `HiC_Contact` | Hi-C contact (cell's own map) | contact | observed/expected contact in the cell's own Hi-C, mapped to K562 5 kb within distance strata | distance-stratified quantile mapping to K562 5 kb | +0.483 |  |

Coefficients are standardized (after the `sign(x)·log1p(|x|)` transform and standardization).
