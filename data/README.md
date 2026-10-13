# Data

## Shipped here (`data/crispr/`)

CRISPR element–gene pairs from the ENCODE CRISPR benchmark used to train and evaluate ENCODE-rE2G
(Gschwind et al.), with CisLink's harmonized input features. ENCODE data are released without use restrictions;
please cite the original studies (below) when using them.

| File | Contents |
|---|---|
| `K562_training.features.tsv.gz` | The 10,250 K562 training pairs (469 positives): label, chromosome fold (1–5) used for tuning, element and promoter coordinates, the 38 model features |
| `{K562,HCT116,GM12878}_heldout.pairs.tsv.gz` | Held-out pairs: label, gene, source study (K562), element–TSS distance, ENCODE-rE2G / ENCODE-rE2G Extended / ABC predictions (where available), distance baseline, and the reference CisLink scores used in the manuscript |
| `{K562,HCT116,GM12878}_heldout.features.tsv.gz` | Feature rows of those pairs, one per element and transcript of the gene; a pair's score is the maximum over its rows |

Held-out K562 pairs come from five studies not used in training (DC-TAP, Klann et al., Morris et al., Reilly et al.,
Xie et al.; 1,895 pairs, 115 positives, after excluding 7 pairs whose gene has no protein-coding TSS). HCT116 has
373 pairs (34 positives) and GM12878 68 pairs (16 positives). For HCT116, ENCODE-rE2G is the better of its two
available predictions.

Coordinates are GRCh38. Features follow `docs/features.md`. Distance baseline = negative element–TSS distance.

## Not redistributed

The manuscript also evaluates CRISPR screens of bone mineral density GWAS variants in hFOB cells: CRISPRi of
89 elements (Conery et al. 2025) and STING-seq of 332 variants (Kahlani et al. 2026), and GWAS gene-prioritization
benchmarks (UK Biobank silver standard, github.com/Deylab999MSKCC/e2g-benchmarking). Their labels are not
redistributed here; obtain them from the original publications.

## Sources

- Gschwind, A. R. *et al.* An encyclopedia of enhancer–gene regulatory interactions in the human genome (ENCODE-rE2G).
- Fulco, C. P. *et al.* (2019); Nasser, J. *et al.* (2021); Morris, J. A. *et al.* (2023); Reilly, S. K. *et al.* (2021);
  Klann, T. S. *et al.* (2021); Xie, S. *et al.* (2019); Ray, J. *et al.* (2025); Guckelberger, P. *et al.* (2024).
