# Data

## Shipped here (`data/crispr/`)

Held-out CRISPR element–gene pairs from the ENCODE CRISPR benchmark (Gschwind et al.), with CisLink's harmonized
input features. ENCODE data are released without use restrictions;
please cite the original studies (below) when using them.

| File | Contents |
|---|---|
| `{K562,HCT116,GM12878,WTC11}_heldout.pairs.tsv.gz` | Held-out pairs: label, gene, source study (K562), element–TSS distance, ENCODE-rE2G / ENCODE-rE2G Extended / ABC predictions (where available), distance baseline, and the reference CisLink v1.3 scores used in the manuscript |
| `{K562,HCT116,GM12878,WTC11}_heldout.features.tsv.gz` | Feature rows of those pairs, one per element and transcript of the gene; a pair's score is the maximum over its rows |

Held-out K562 pairs come from five studies not used in training (DC-TAP, Klann et al., Morris et al., Reilly et al.,
Xie et al.; 1,895 pairs, 115 positives, after excluding 7 pairs whose gene has no protein-coding TSS). HCT116 has
373 pairs (34 positives), GM12878 68 pairs (16 positives) and WTC11 DC-TAP 1,908 pairs (15 positives). For HCT116,
ENCODE-rE2G is the better of its two available predictions. WTC11 has no histone marks, so only CisLink-V2G is scored
there (CisLink-E2G columns are empty).

Inputs follow the manuscript's rules: accessibility as a MACS2 -log10 P track (K562 and GM12878: Roadmap DNase-seq;
HCT116 and WTC11: ENCODE DNase-seq alignments, biological replicates pooled, processed by the Roadmap procedure,
`cislink pvalue --assay dnase`); expression from poly(A) RNA-seq, the mean over biological replicates (ENCODE RSEM;
WTC11 the mean of eight GSE314132 libraries quantified with the ENCODE pipeline); Hi-C as each cell's 5-kb O/E. The
v1.3 cross-epigenome correlations (`F1_dnase`, `F1_h3k27ac`) are the per-element values used in the manuscript
(computed from the Roadmap bigWig tracks); recomputing them from the released 25-bp reference (`cislink add-f1`) changes
them slightly (Pearson r >= 0.9995) and the AUPRCs by at most 0.0023. `DxA` is computed by the package.

Coordinates: K562 and GM12878 keep the hg19 coordinates of the training data; HCT116 and WTC11 are GRCh38. Features
follow `docs/features.md`. Distance baseline = negative element–TSS distance.

## GWAS silver standard (`data/gwas/silver_standard.tsv.gz`)

Every protein-coding candidate gene within 1 Mb of the 568 noncoding UK Biobank credible sets for 31 traits whose one gene within
2 Mb carries an independent coding variant for the same trait (the silver gene; Weeks et al. 2023, Gschwind et al.;
github.com/Deylab999MSKCC/e2g-benchmarking), as in Supplementary Table 6. Each trait is scored in its matched cell type (osteoblast, K562, CD4+ T cell, CD14+
monocyte, HepG2, kidney or IMR-90; Supplementary Table 1 lists the inputs). Columns: trait, cell type scored, credible set,
fine-mapping region, candidate gene, silver-gene flag, CisLink-V2G and CisLink-E2G (maximum over the set's variants
with PIP > 0.1), the nearest-gene baseline (negative distance from the highest-probability variant to the gene's nearest
protein-coding TSS, GENCODE v29 on GRCh38), the published
methods' scores (ENCODE-rE2G, ABC, cS2G, PoPS, L2G; FLAMES on eight traits), and the leave-out flags (exome-burden
silver gene, L2G training gene).

## Not redistributed

The manuscript also evaluates CRISPR screens of bone mineral density GWAS variants in hFOB cells: CRISPRi of
89 elements (Conery et al. 2025) and STING-seq of 332 variants (Kahlani et al. 2026). Their labels are not redistributed
here; obtain them from the original publications.

## Sources

- Gschwind, A. R. *et al.* An encyclopedia of enhancer–gene regulatory interactions in the human genome (ENCODE-rE2G).
- Weeks, E. M. *et al.* (2023) Leveraging polygenic enrichments of gene features to predict genes underlying complex traits and diseases.
- Fulco, C. P. *et al.* (2019); Nasser, J. *et al.* (2021); Morris, J. A. *et al.* (2023); Reilly, S. K. *et al.* (2021);
  Klann, T. S. *et al.* (2021); Xie, S. *et al.* (2019); Ray, J. *et al.* (2025); Guckelberger, P. *et al.* (2024).
