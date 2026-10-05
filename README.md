# rsemflow

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemflow/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Downstream bulk RNA-seq analysis of RSEM `*.genes.results` files, driven by an
ordinary sample metadata table and run from the command line or from R.

```text
*.genes.results + metadata.tsv
        │  rsemflow data import
        ▼
     study/  ──► summaries · annotation · normalization · PCA · correlation
                 DESeq2 · preranked GSEA · GSVA / ssGSEA
        ▼
   TSV tables
```

rsemflow connects established Bioconductor methods (tximport, DESeq2, fgsea,
GSVA, msigdbr) behind one consistent interface. It writes tables only — no
figures, HTML reports, or Shiny apps — so results can be plotted or analysed
further in R, Python, Excel, or Prism.

Two principles guide the design:

- **You state the experimental design.** The model formula, the term to test,
  and the reference levels are always given explicitly. rsemflow never guesses
  a design, removes samples, or interprets results.
- **Stop rather than return a meaningless result.** Mismatched sample IDs,
  inconsistent gene tables, missing values in model variables, and
  rank-deficient (confounded) designs all stop the run with a message.

## Installation

rsemflow needs R ≥ 4.3 and several Bioconductor packages.

### From GitHub (R)

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemflow")
```

`BiocManager` resolves the Bioconductor dependencies. Human and mouse data are
supported.

### From a source checkout (with the command-line launcher)

On Linux, macOS, WSL, or an HPC login node:

```bash
git clone https://github.com/DaikiKumakura/rsemflow.git
cd rsemflow
bash install_local.sh
```

The installer installs missing dependencies (`scripts/install_dependencies.R`),
runs `R CMD INSTALL .`, links the `rsemflow` launcher into
`${RSEMFLOW_BIN:-$HOME/.local/bin}`, and prints the installed version. Set
`RSEMFLOW_SKIP_DEPS=1` to skip the dependency step. Add the launcher directory
to `PATH` if needed:

```bash
export PATH="$HOME/.local/bin:$PATH"
rsemflow --help
```

On Windows without WSL, call the same CLI through `Rscript`:

```bash
Rscript -e "rsemflow::cli_main()" data import --i-rsem-dir data --m-metadata metadata.tsv --o-study study
```

## Input

```text
project/
├── data/
│   ├── S01.genes.results
│   ├── S02.genes.results
│   └── ...
└── metadata.tsv
```

**RSEM files.** The sample ID is the file name without `.genes.results`. Each
file needs the standard columns `gene_id`, `transcript_id(s)`, `length`,
`effective_length`, `expected_count`, `TPM`, and `FPKM`. All files must contain
the same genes in the same order, without duplicated `gene_id`.

**Metadata.** A plain TSV or CSV (by file extension) with a `sample_id` column;
every other column is optional.

```text
sample_id	condition	batch	age
S01	Control	A	51
S02	Control	B	43
S03	Drug	A	48
S04	Drug	B	55
```

Every RSEM sample must appear in the metadata and vice versa. `sample_id` is
always read as text, so `001` stays `001`. A column whose values all parse as
numbers becomes numeric; any other column becomes categorical. Override the
guess with `--p-categorical COLUMN` (e.g. time points `0, 6, 24` as three
groups) or `--p-numeric COLUMN`.

## Command-line interface

```text
rsemflow <module> <action> [options]
```

| Module | Actions |
| --- | --- |
| `data` | `import`, `summarize` |
| `annotation` | `ensembl`, `gtf` |
| `expression` | `normalize`, `pca`, `correlation` |
| `differential` | `deseq2` |
| `enrichment` | `gsea`, `gsva` |

Option prefixes tell you what each option is: `--i-*` inputs, `--m-*`
metadata, `--p-*` parameters, `--o-*` outputs. Every level has help:
`rsemflow --help`, `rsemflow differential --help`,
`rsemflow differential deseq2 --help`.

### A typical workflow

```bash
# 1. Import and validate
rsemflow data import --i-rsem-dir data/ --m-metadata metadata.tsv --o-study study/

# 2. Expression-level summaries (not alignment QC)
rsemflow data summarize --i-study study/ --o-summary results/summary/

# 3. Sample structure
rsemflow expression pca --i-study study/ --o-pca results/pca/
rsemflow expression correlation --i-study study/ --o-correlation results/correlation.tsv

# 4. Differential expression
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ batch + condition' \
  --p-term condition \
  --p-reference condition::Control \
  --o-differential results/differential/

# 5. Pathways
rsemflow enrichment gsea \
  --i-differential results/differential/condition_Drug_vs_Control.tsv \
  --p-species human --p-collection hallmark \
  --o-enrichment results/gsea.tsv
rsemflow enrichment gsva \
  --i-study study/ --p-species human --p-collection hallmark \
  --o-scores results/gsva.tsv

# 6. Gene annotation (prefer the GTF used to build the RSEM reference)
rsemflow annotation gtf --i-study study/ --i-gtf genes.gtf.gz --o-annotation results/annotation.tsv
```

A complete run on the bundled synthetic data (6 samples, 200 genes; no
downloads needed):

```bash
bash inst/extdata/example/run_example.sh example-output
```

### Differential-expression designs

`--p-formula` is the full model, `--p-term` chooses which results to write, and
`--p-reference COLUMN::LEVEL` (repeatable) sets factor reference levels.

| Design | Options | Tables written |
| --- | --- | --- |
| Two groups | `--p-formula '~ condition' --p-term condition --p-reference condition::Control` | `condition_Drug_vs_Control.tsv` |
| Batch-adjusted | `--p-formula '~ batch + condition' --p-term condition ...` | as above |
| Paired | `--p-formula '~ patient + condition' --p-term condition --p-reference condition::Pre` | `condition_Post_vs_Pre.tsv` |
| Several groups | `--p-term treatment --p-reference treatment::Vehicle` | every level vs `Vehicle`; add `--p-contrast treatment::DrugB::DrugA` for others |
| Numeric covariate | `--p-formula '~ sex + age + score' --p-term score` | `score.tsv` (per-unit change) |
| Interaction | `--p-formula '~ genotype * treatment' --p-term genotype:treatment` | one table per interaction coefficient |
| Simple effect | interaction formula + `--p-term treatment --p-contrast treatment::Drug::Vehicle --p-at genotype::KO` | `treatment_Drug_vs_Vehicle_at_genotype_KO.tsv` |
| Omnibus (LRT) | `--p-test lrt --p-formula '~ patient + time' --p-reduced '~ patient' --p-term time` | `time_omnibus.tsv` |

Notes:

- When the term interacts with another factor (`~ genotype * treatment`), the
  plain `--p-term treatment` comparisons are effects at the reference level of
  the other factor. rsemflow warns about this; use `--p-at` for other levels.
- `--p-at` takes exactly one `--p-at` and one `--p-contrast`, and the contrast
  denominator must be the reference level.
- LRT tables have `log2FoldChange` and `lfcSE` set to `NA`: the omnibus test
  has no single direction. They cannot be used for preranked GSEA.
- Log fold changes are the ordinary DESeq2 estimates (no shrinkage).
- The design matrix must be full rank. rsemflow does not judge whether you
  have enough biological replicates.

The output directory also contains `design-matrix.tsv` (the model matrix
used), `contrasts.tsv` (one row per table, with the DESeq2 coefficient), and
`analysis-info.tsv` (formula, term, test, reduced formula, alpha).

### Gene sets

`--p-collection` accepts `hallmark`, `reactome`, `go-bp`, `go-mf`, or `go-cc`,
or a raw MSigDB code with `--p-subcollection` (for example
`--p-collection C2 --p-subcollection CP:KEGG_MEDICUS`). Gene sets come from
`msigdbr`. Human uses the human MSigDB; mouse uses the mouse MSigDB, where the
aliases map to `MH`, `M2`, and `M5` (raw codes must be mouse codes). Use `--i-genesets file.gmt` for custom
sets.

Genes are matched by stable Ensembl gene ID (`ENSG00000141510.18` →
`ENSG00000141510`). For GSEA, when several rows share a stable ID the one with
the largest absolute statistic is kept; for GSVA their expression is averaged.
`fgseaMultilevel` estimates p-values by sampling, so GSEA uses a fixed seed
(`--p-seed`, default 42) and repeated runs give identical tables.

## Outputs

| Command | Files |
| --- | --- |
| `data import` | `study.rds`, `metadata.tsv`, `samples.tsv`, `expected_count.tsv`, `tpm.tsv`, `effective_length.tsv`, `study-info.tsv` |
| `data summarize` | `samples.tsv` (total counts, detected genes, TPM quantiles), `genes.tsv` (mean/median count and TPM, detection) |
| `annotation ensembl` / `gtf` | one row per RSEM gene; the original `gene_id` is kept alongside `ensembl_gene_id`, symbol, and status |
| `expression normalize` | gene × sample matrix (`vst`, `normalized-counts`, or `log2-tpm`) |
| `expression pca` | `scores.tsv`, `loadings.tsv`, `variance.tsv`, `analysis-info.tsv` (VST, top 500 variable genes, centred, unscaled by default) |
| `expression correlation` | sample × sample matrix (VST + Pearson by default) |
| `differential deseq2` | one table per comparison plus the three manifest files above |
| `enrichment gsea` | fgsea columns, `leadingEdge` as `GENE1;GENE2`, plus collection and MSigDB version; sorted by `padj` then `|NES|` |
| `enrichment gsva` | pathway × sample scores from `log2(TPM + 1)` |

RSEM sets `effective_length` to 0 for genes shorter than the fragment length.
Those genes have zero counts; rsemflow keeps the original values in the study
and replaces the zeros with 1 only when building the DESeq2 object, as the
DESeq2 vignette recommends.

## Using rsemflow from R

Every CLI action calls an exported R function:

```r
library(rsemflow)

study <- read_rsem_study("data", "metadata.tsv")
write_study(study, "study")

summary <- summarize_study(study)
pca <- compute_pca(study)
cor_mat <- compute_correlation(study)

de <- differential_deseq2(
  study,
  formula = ~ batch + condition,
  term = "condition",
  reference = "condition::Control"
)
gsea <- run_gsea(de$results$condition_Drug_vs_Control, species = "human")
gsva <- run_gsva(study, species = "human")
```

## Scope

Included: RSEM gene-level import, expression summaries, annotation,
normalization, PCA, correlation, DESeq2, GSEA, GSVA/ssGSEA.

Not included: plotting, reports, automatic sample removal or design selection,
LFC shrinkage, alignment-level QC, transcript-level or splicing analysis,
testing of GSVA scores between groups.

## Development

```bash
make deps      # install dependencies
make document  # regenerate man/ and NAMESPACE (roxygen2)
make test      # testthat suite
make check     # R CMD build + R CMD check
```

## Citation

rsemflow only coordinates other methods; please cite those you use
(`citation("rsemflow")` lists RSEM, tximport, DESeq2, fgsea, GSVA, and MSigDB).

## License

MIT © Daiki Kumakura
