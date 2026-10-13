# Data

## Shipped here (`data/crispr/`)

Held-out CRISPR element–gene pairs from the ENCODE CRISPR benchmark (Gschwind et al.), with CisLink's harmonized
input features. ENCODE data are released without use restrictions;
please cite the original studies (below) when using them.

| File | Contents |
|---|---|
| `{K562,HCT116,GM12878}_heldout.pairs.tsv.gz` | Held-out pairs: label, gene, source study (K562), element–TSS distance, ENCODE-rE2G / ENCODE-rE2G Extended / ABC predictions (where available), distance baseline, and the reference CisLink scores used in the manuscript |
| `{K562,HCT116,GM12878}_heldout.features.tsv.gz` | Feature rows of those pairs, one per element and transcript of the gene; a pair's score is the maximum over its rows |

Held-out K562 pairs come from five studies not used in training (DC-TAP, Klann et al., Morris et al., Reilly et al.,
Xie et al.; 1,895 pairs, 115 positives, after excluding 7 pairs whose gene has no protein-coding TSS). HCT116 has
373 pairs (34 positives) and GM12878 68 pairs (16 positives). For HCT116, ENCODE-rE2G is the better of its two
available predictions.

Coordinates are GRCh38. Features follow `docs/features.md`. Distance baseline = negative element–TSS distance.

## GWAS silver standard (`data/gwas/silver_standard.tsv.gz`)

Every candidate gene within 1 Mb of the 568 noncoding UK Biobank credible sets for 31 traits whose one gene within
2 Mb carries an independent coding variant for the same trait (the silver gene; Weeks et al. 2023, Gschwind et al.;
github.com/Deylab999MSKCC/e2g-benchmarking), as in Supplementary Table 6: trait, cell type scored, credible set,
fine-mapping region, candidate gene, silver-gene flag, CisLink-V2G and CisLink-E2G (maximum over the set's variants
with PIP > 0.1), the nearest-gene baseline (negative TSS distance to the highest-probability variant), the published
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
