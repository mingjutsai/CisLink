# Computing features for a new cell type

`cislink build` turns a cell type's raw data into the harmonized feature table the models score (CisLink v1.3:
CisLink-E2G 40 features, CisLink-V2G 16). It reproduces the pipeline used for the manuscript: on the eBMD variants of
chromosome 22 (13,243 variant–transcript rows, osteoblasts), the cell-derived features matched the published matrix
(relative difference < 1e-4). The v1.3 cross-epigenome correlations (F1) are computed by `cislink.roadmap` from the
Roadmap activity reference; on the manuscript's evaluation sets this code reproduces the reference F1 values exactly
(see "Cross-epigenome correlations" below).

## Install

```bash
pixi install -e pipeline        # adds bedtools, pyBigWig, hic-straw, PyYAML (Linux, macOS)
pixi install -e macs2           # CrossMap (lifts elements to hg19 for the F1 features) and MACS2
pixi run cislink reference download --dest /data/cislink_ref    # Roadmap activity reference, 16.3 GB (once)
pixi run -e pipeline cislink build --help
```

## Accessibility from a BAM

The models read accessibility as a MACS2 -log10 P track, the format of their training data (Roadmap consolidated
DNase-seq). Use the published -log10 P track where one exists (Roadmap DNase-seq, ENCODE ATAC-seq "signal p-value");
do not use read-depth or fold-change tracks. When a cell has only aligned reads, `cislink pvalue` makes the track from
the BAM, with one procedure per assay:

- **DNase-seq: the Roadmap procedure** that produced the training track (reconstructed; on K562 it reproduces the
  Roadmap E123 track, element-level Spearman 0.998). Tags: MAPQ ≥ 30, primary, not duplicate-flagged, read 1 of
  pairs, truncated to their 5' 36 bp, at most 50 million (fixed-seed subsample with GNU `shuf`; on macOS install
  coreutils and put its `gnubin` directory on the PATH). Control: Roadmap's simulated uniform
  DNase control lifted to GRCh38. Biological replicates are pooled: pass every replicate's BAM, and their tags are combined
  before the 50-million cap, as in Roadmap's consolidated tracks. `macs2 callpeak -c control -f BED -g hs -p 0.01 --nomodel --extsize 146 -B --SPMR`
  (equal to Roadmap's MACS2 2.0.10 `--shiftsize 73`), then `macs2 bdgcmp -m ppois -S <min(tags, control) / 1e6>`.
- **ATAC-seq: the ENCODE ATAC-seq pipeline.** Filtered, deduplicated reads; Tn5 insertion sites (+4/−5);
  `macs2 callpeak -f BED --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01` without a control; and
  `macs2 bdgcmp -m ppois -S <cut sites / 1e6>`. Duplicates are marked first when the BAM header shows no marking
  step (`--mark-duplicates auto`). Several BAMs (replicates) are deduplicated one by one and their cut sites pooled, as
  the ENCODE pipeline's pooled-replicate track.

```bash
pixi install -e macs2           # MACS2 2.2.7.1 and CrossMap (Linux, macOS x86-64)
# DNase-seq: fetch and lift the Roadmap control once (about 1.6 GB), then process each BAM
pixi run -e pipeline cislink roadmap-control --crossmap .pixi/envs/macs2/bin/CrossMap -o roadmap_dnase_control.hg38.bed
pixi run -e pipeline cislink pvalue rep1.bam rep2.bam --assay dnase --control roadmap_dnase_control.hg38.bed \
    --macs2 .pixi/envs/macs2/bin/macs2 -o sample_DNase_pval.bedGraph
# ATAC-seq
pixi run -e pipeline cislink pvalue sample.nodup.bam --assay atac --macs2 .pixi/envs/macs2/bin/macs2 -o sample_ATAC_pval.bedGraph
```

Validation: the DNase procedure run on the HCT116 DNase-seq alignments (ENCFF304MDQ) reproduces the track used in the
manuscript (chr1 and chr21: 255 of 296 million bases differ, by at most 0.14; Spearman 1.0000 over 500-bp windows),
and `cislink roadmap-control` reproduces the manuscript's lifted control byte for byte; on osteoblast ATAC-seq the ATAC procedure agrees with the ENCODE pipeline's own
p-value track (Spearman 0.996 over the candidate elements), and on myotube ATAC-seq it reproduces the manuscript's
track (element Spearman 0.99999). Applying the DNase procedure to ATAC-seq gave the same results as the ATAC procedure,
so each assay keeps its own procedure. Use the track as `accessibility` and as `accessibility_pvalue` in the
configuration.

## Expression

Use poly(A) RNA-seq quantified with the ENCODE long-RNA pipeline (STAR, RSEM; GENCODE v29), the format of the
training data (K562 poly(A) RNA-seq, mean TPM of the two biological replicates). Total RNA-seq is not equivalent:
intronic and pre-mRNA reads shift the TPM of long and short genes in gene-specific ways that rescaling cannot undo. In
the manuscript every cell used the mean TPM over its biological replicates or eligible samples (ENCODE RSEM files, or
RNA-seq quantified from public reads with the library strand inferred per sample, e.g. from STAR
`--quantMode GeneCounts`, and adapters trimmed where present). Give the averaged table as the quantification. Check a
self-quantified sample against a published quantification of the same library before use.

## Inputs (GRCh38)

| Input | Used for | Format |
|---|---|---|
| GENCODE annotation (v29 in the manuscript) | candidate transcripts and genes, gene context | GTF (`.gtf` or `.gtf.gz`) |
| Chromatin tracks: `accessibility` (ATAC-seq or DNase-seq), 11 histone marks, `CTCF` | `{mark}_Enh`, `{mark}_Prom` | bedGraph (through `bedtools`) or bigWig; signal or -log10 P, as the cell's tracks come |
| Candidate elements of the cell | quantile mapping of chromatin to K562 | BED (chrom, start, end, activity), or an accessibility -log10 P bedGraph to call them |
| Transcript and gene quantifications | `Tx_expression`, `Gene_expression` | ENCODE RSEM TSVs (`transcript_id` / `gene_id`, `TPM`) |
| The cell's Hi-C (CisLink-E2G only) | `HiC_Contact` | `.hic` |
| Multi-tissue Hi-C megamap (CisLink-V2G) | `HiC_megamap_oe` | `.hic` path or URL (default: the ENCODE-rE2G megamap) |
| Roadmap activity reference (both models) | `F1_dnase`, `F1_h3k27ac` | directory from `cislink reference download` (`--roadmap-ref`, `$CISLINK_ROADMAP_REF` or `roadmap_reference:`) |
| CrossMap (both models) | GRCh38 → hg19 lift of elements for F1 | executable (`pixi install -e macs2`; `--crossmap` or `crossmap:`) |

CisLink-V2G needs only accessibility, expression, the megamap and the Roadmap reference; leave the histone tracks and
`hic` out to build V2G features only (E2G is then not scored). The reference is the same for every cell type.

## Input formats

| Input | Accepted |
|---|---|
| Variants (`--variants`) | TSV with columns `SNPID`, `chr`, `pos` (1-based), or without a header: chromosome, position and optionally an ID (default `chr:pos`) |
| Elements (`--elements`) | BED: chrom, start, end, name (0-based, half-open) |
| Expression | ENCODE RSEM output (`transcript_id` / `gene_id` with a `TPM` column), or any table with the Ensembl ID in the first column and the value in the last (header optional). Version suffixes are ignored (`ENSG…​.11` = `ENSG…`); PAR_Y rows are dropped; IDs absent from the table get 0 |
| Chromatin tracks | bedGraph (chrom, start, end, value; `track`/`browser`/`#` lines skipped; `.gz` allowed) or bigWig (`.bw`, `.bigWig`) |
| Candidate elements | BED: chrom, start, end, activity |

**Chromosome naming is detected per file.** Variants, elements, tracks and Hi-C maps may use `chr1` or `1`; names are
converted to the annotation's convention internally. This matters because a bedGraph named `1` scored against
`chr1` regions would otherwise give all-zero signal without an error. bedGraph tracks must agree with each other.

`cislink validate` runs the same preflight that `cislink build` runs first: required keys, readable files, track
names (an `accessibility` track is required; `DNase` and `ATAC` are accepted as its name), the first data row of each bedGraph, chromosome naming, bedtools, the
stable Ensembl IDs shared by the annotation and the expression tables, duplicated variant IDs, and element
coordinates. Problems are reported before any long step starts.

```bash
pixi run -e pipeline cislink validate --variants variants.tsv --config osteoblast.yaml
```

## Configuration

```yaml
# osteoblast.yaml
gtf: gencode.v29.annotation.gtf
transcript_quantification: osteoblast_RNAseq_transcripts.tsv
gene_quantification: osteoblast_RNAseq_genes.tsv
tracks:
  accessibility: osteoblast_ATAC_pval.bedGraph   # ATAC-seq or DNase-seq
  H3K27ac: osteoblast_H3K27ac.bigWig
  H3K4me1: osteoblast_H3K4me1.bigWig
  # ... H4K20me1, H3K9me3, H3K9ac, H3K79me2, H3K4me3, H3K4me2, H3K36me3, H3K27me3, H2A.Z, CTCF
candidates: osteoblast_candidates.bed       # or: accessibility_pvalue: osteoblast_ATAC_pval.bedGraph
hic: osteoblast.hic                         # omit for V2G only
hic_resolution: 2500                        # the map's finest reliable resolution (5000 for most maps)
hic_norm: KR
# map_hic: true                             # default: true for maps finer than 5 kb, false at 5 kb
megamap: tissues.hic                        # local copy speeds things up; default is the public URL
roadmap_reference: /data/cislink_ref/roadmap_activity_v1   # or --roadmap-ref / $CISLINK_ROADMAP_REF
crossmap: .pixi/envs/macs2/bin/CrossMap
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
# check the inputs first (build runs the same checks)
pixi run -e pipeline cislink validate --variants variants.tsv --config osteoblast.yaml
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
4. **Hi-C.** Observed/expected contact between the element-midpoint bin and the TSS bin. Maps at the training resolution (5 kb)
   are used directly, with the map's own normalization (KR; SCALE for intact Hi-C, an equivalent matrix balancing).
   Finer maps (`hic_resolution` below 5000) keep their resolution and are mapped onto K562 5-kb O/E within 16 distance
   strata, using the cell's own distribution sampled from 20 random 2-Mb windows per autosome (`map_hic`, on by
   default only for such maps). The megamap is read at 5 kb (SCALE normalization). Missing contact is scored as 0.
5. **Distance and gene context.** Distance (element midpoint to TSS), log distance, element length, power-law contact
   ((max(d, 5 kb) / 5 kb)^-1.024), protein-coding genes strictly between element and TSS, the gene's distance rank
   among protein-coding genes (1 = nearest), nearest-gene status, and ABC's list of ubiquitously expressed genes.
   `DxA` = log10 distance × harmonized element accessibility (computed when scoring).
6. **Cross-epigenome correlations (F1).** See below.

The K562 references (`src/cislink/reference/`) come from ENCODE K562 data processed the same way.

## Cross-epigenome correlations (F1)

`F1_dnase` and `F1_h3k27ac` are the Spearman correlations, across Roadmap reference epigenomes, between the element's
activity and the gene's expression: the element's Roadmap consolidated MACS2 −log10 P signal (DNase-seq, 28
epigenomes; H3K27ac, 47 epigenomes) and log(1 + RPKM) of the gene in the Roadmap protein-coding RNA-seq matrix. The
epigenomes are every one with RNA-seq and an observed (not imputed) track of the mark, minus E000 (pooled reference
RNA) and E123 (K562, the training cell). Imputed tracks are not used because they are predicted from the same
epigenome's other marks, which would make the correlation partly circular.

- **Reference.** `cislink reference download --dest DIR` fetches the Roadmap activity reference (Zenodo; 16.3 GB):
  the 75 tracks binned at 25 bp on hg19, the RNA matrix and the UCSC hg38ToHg19 chain, and verifies every file (md5).
  `cislink reference verify --dest DIR/roadmap_activity_v1` re-checks it.
- **Computation** (`cislink build`, or `cislink add-f1` for an existing feature table with `Enh` and `TX`
  columns). Elements are lifted from GRCh38 to hg19 with CrossMap (one-to-one only); the element signal is the sum of
  its 25-bp bins, edge bins prorated by the covered fraction; the gene is identified through the transcript's GTF gene
  ID (symbol as fallback). A pair gets 0 when the element does not lift, the gene has no Roadmap RNA-seq, or a vector
  is constant.
- **Same for every cell type.** F1 describes the element–gene pair across tissues; cell-specific evidence comes from the
  cell's own accessibility, expression and Hi-C. The Roadmap epigenomes with RNA-seq include no osteoblasts or
  myotubes.
- **Agreement with the manuscript.** The manuscript computed F1 from the bigWig tracks per element. The 25-bp
  reference reproduces those values with Pearson r ≥ 0.9995 in every evaluation set; with the reference values the
  AUPRC of the models changes by at most 0.0023 and top-1 accuracy is unchanged. `cislink.roadmap` reproduces the
  reference values of all 18 sets (max |difference| 5e-7, the rounding of the stored values).

## Reading the output

Rank genes per variant or element with `cislink rank`. Scores are probabilities on the scale of the K562 CRISPR
training pairs; compare them with the recall-based thresholds (E2G 0.195, V2G 0.147 at 70% recall) rather than
reading them as absolute probabilities (notebook 04).
