# Computing features for a new cell type

`cislink build` turns a cell type's raw data into the harmonized feature table the models score. It reproduces
the pipeline used for the manuscript: on the eBMD variants of chromosome 22 (13,243 variant–transcript rows,
osteoblasts), every one of the 38 features and both models' scores match the published matrix (relative
difference < 1e-4; scores < 1e-6).

## Install

```bash
pixi install -e pipeline        # adds bedtools, pyBigWig, hic-straw, PyYAML (Linux, macOS)
pixi run -e pipeline cislink build --help
```

## Inputs (GRCh38)

| Input | Used for | Format |
|---|---|---|
| GENCODE annotation (v29 in the manuscript) | candidate transcripts and genes, gene context | GTF (`.gtf` or `.gtf.gz`) |
| Chromatin tracks: accessibility (`DNase`, ATAC-seq or DNase-seq), 11 histone marks, `CTCF` | `{mark}_Enh`, `{mark}_Prom` | bedGraph (through `bedtools`) or bigWig; signal or -log10 P, as the cell's tracks come |
| Candidate elements of the cell | quantile mapping of chromatin to K562 | BED (chrom, start, end, activity), or an accessibility -log10 P bedGraph to call them |
| Transcript and gene quantifications | `Tx_expression`, `Gene_expression` | ENCODE RSEM TSVs (`transcript_id` / `gene_id`, `TPM`) |
| The cell's Hi-C (CisLink-E2G only) | `HiC_Contact` | `.hic` |
| Multi-tissue Hi-C megamap (CisLink-V2G) | `HiC_megamap_oe` | `.hic` path or URL (default: the ENCODE-rE2G megamap) |

CisLink-V2G needs only accessibility, expression and the megamap; leave the histone tracks and `hic` out to build
V2G features only (E2G is then not scored).

## Configuration

```yaml
# osteoblast.yaml
gtf: gencode.v29.annotation.gtf
transcript_quantification: osteoblast_RNAseq_transcripts.tsv
gene_quantification: osteoblast_RNAseq_genes.tsv
tracks:
  DNase: osteoblast_ATAC_pval.bedGraph      # accessibility, ATAC-seq or DNase-seq
  H3K27ac: osteoblast_H3K27ac.bigWig
  H3K4me1: osteoblast_H3K4me1.bigWig
  # ... H4K20me1, H3K9me3, H3K9ac, H3K79me2, H3K4me3, H3K4me2, H3K36me3, H3K27me3, H2A.Z, CTCF
candidates: osteoblast_candidates.bed       # or: accessibility_pvalue: osteoblast_ATAC_pval.bedGraph
hic: osteoblast.hic                         # omit for V2G only
hic_resolution: 2500                        # the map's finest reliable resolution (5000 for most maps)
hic_norm: KR
map_hic: true                               # map O/E onto K562 within distance strata
megamap: tissues.hic                        # local copy speeds things up; default is the public URL
bedtools: bedtools
threads: 8
```

See `examples/config.example.yaml`.

The first run computes the cell's chromatin quantiles (over its candidate elements) and, for E2G, its Hi-C quantiles,
and saves them next to the output (`OUT.chromatin_quantiles.tsv`, `OUT.hic_quantiles.tsv`). Point later runs of the
same cell at them to skip that step:

```yaml
chromatin_quantiles: features.tsv.chromatin_quantiles.tsv
hic_quantiles: features.tsv.hic_quantiles.tsv
```

Quantiles are kept at 6 significant digits, as in the manuscript pipeline, so results reproduce exactly.

## Run

```bash
# GWAS variants (TSV: SNPID, chr, pos) -> every protein-coding transcript TSS within 1 Mb
pixi run -e pipeline cislink build --variants variants.tsv --config osteoblast.yaml --score -o features.tsv
# candidate elements (BED: chrom, start, end, name)
pixi run -e pipeline cislink build --elements elements.bed --config osteoblast.yaml --score -o features.tsv
# one row per variant (or element) and gene: maximum over transcripts, ranks, threshold flags
pixi run cislink rank features.tsv -o genes.tsv
```

## What each step does

1. **Pairs.** A variant's element is the variant ±250 bp; an element is used as given. The promoter is each
   protein-coding transcript TSS ±250 bp (a transcript whose TSS equals that of the preceding transcript in the GTF
   is skipped). One row per transcript; a pair's score is the maximum over its rows.
2. **Chromatin.** sum(signal × overlapping bp) over the element and the promoter. The cell's distribution of per-bp
   mean signal over its genome-wide candidate elements (1,001 quantiles) is mapped onto K562's; each sum becomes a
   per-bp mean, is mapped, and is multiplied by the region length again. Candidates called from an accessibility
   -log10 P bedGraph: bins at or above a threshold merged across gaps ≤ 100 bp, with the threshold giving closest to
   150,000 regions.
3. **Expression.** Transcript TPM matched on the stable Ensembl ID; gene TPM through the gene symbol's GTF gene ID;
   both rescaled to protein-coding TPM (1e6 / the protein-coding TPM sum).
4. **Hi-C.** Observed/expected contact between the element-midpoint bin and the TSS bin. With `map_hic`, the cell's
   O/E is mapped onto K562 5-kb O/E within 16 distance strata, using the cell's own distribution sampled from 20
   random 2-Mb windows per autosome. The megamap is read at 5 kb (SCALE normalization). Missing contact is scored as 0.
5. **Distance and gene context.** Distance (element midpoint to TSS), log distance, element length, power-law contact
   ((max(d, 5 kb) / 5 kb)^-1.024), protein-coding genes strictly between element and TSS, the gene's distance rank
   among protein-coding genes (1 = nearest), nearest-gene status, and ABC's list of ubiquitously expressed genes.

The K562 references (`src/cislink/reference/`) come from ENCODE K562 data processed the same way.

## Reading the output

Rank genes per variant or element with `cislink rank`. Scores are probabilities on the scale of the K562 CRISPR
training pairs; compare them with the recall-based thresholds (E2G 0.180, V2G 0.139 at 70% recall) rather than
reading them as absolute probabilities (notebook 04).
