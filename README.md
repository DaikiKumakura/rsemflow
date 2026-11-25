# rsemflow

<!-- badges: start -->
[![R-universe](https://daikikumakura.r-universe.dev/badges/rsemflow)](https://daikikumakura.r-universe.dev/rsemflow)
[![R-CMD-check](https://github.com/DaikiKumakura/rsemflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemflow/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210477.svg)](https://doi.org/10.5281/zenodo.24210477)
<!-- badges: end -->

**Standard bulk RNA-seq results from RSEM files and a sample table, as TSV
tables, from the command line or from R.**

## Installation

```r
install.packages("BiocManager")
install.packages("rsemflow", repos = c(
  "https://daikikumakura.r-universe.dev",
  BiocManager::repositories()
))
```

Or from GitHub: `BiocManager::install("DaikiKumakura/rsemflow")`. rsemflow
needs R 4.3 or later. To get the `rsemflow` command on Linux, macOS, WSL or an
HPC login node, clone the repository and run `bash install_local.sh`; on
Windows, call the same commands as `Rscript -e "rsemflow::cli_main()" data
import ...`. See [Installation](https://daikikumakura.github.io/rsemflow/articles/installation.html).

## Quick start

The package includes a small synthetic study (6 samples, 200 genes), so this
runs as it is, with no downloads:

```r
library(rsemflow)
ex <- system.file("extdata", "example", package = "rsemflow")

study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
summary <- summarize_study(study)
pca <- compute_pca(study)
de <- differential_deseq2(
  study,
  formula = ~ batch + condition,
  term = "condition",
  reference = "condition::Control"
)
head(de$results$condition_Drug_vs_Control)
```

The same from the command line, writing every table to a folder:

```bash
bash inst/extdata/example/run_example.sh example-output
```

## Documentation

[Function reference](https://daikikumakura.github.io/rsemflow/reference/) · [Changelog](https://daikikumakura.github.io/rsemflow/news/) · [Releases](https://github.com/DaikiKumakura/rsemflow/releases)

## Workflow

<details>
<summary>Show details</summary>

rsemflow takes RSEM `*.genes.results` files and an ordinary sample metadata
table and produces the usual downstream results: expression summaries,
normalized expression, PCA, sample correlation, DESeq2 differential
expression, preranked GSEA, GSVA/ssGSEA scores and gene annotation. It
connects established Bioconductor methods (tximport, DESeq2, fgsea, GSVA,
msigdbr) behind one interface and writes tables only, so the results can be
plotted or analysed further in R, Python, Excel or Prism.

- **You state the design.** The model formula, the term to test and the
  reference levels are always given explicitly; rsemflow never guesses a
  design, removes samples or interprets results.
- **It stops rather than return a meaningless result.** Mismatched sample IDs,
  inconsistent gene tables, missing model values and confounded designs stop
  the run with a message.
- **Human and mouse,** with MSigDB gene sets for both.
- **Checked on real data:** on ENCODE knockdown experiments, each knockdown
  target comes out among the most significant decreases in its own knockdown.

**Documentation:** <https://daikikumakura.github.io/rsemflow/> ·
[Get started](https://daikikumakura.github.io/rsemflow/articles/rsemflow.html) ·
[Tutorial](https://daikikumakura.github.io/rsemflow/articles/tutorial.html) ·
[Designs](https://daikikumakura.github.io/rsemflow/articles/designs.html) ·
[Installation](https://daikikumakura.github.io/rsemflow/articles/installation.html)

```text
  data/*.genes.results + metadata.tsv
                 │
   data import                 validate and store the study
                 │
   data summarize              counts, detected genes, TPM quantiles
   expression pca / correlation  sample structure
   expression normalize        VST, normalized counts, log2 TPM
   differential deseq2         one table per comparison
   enrichment gsea / gsva      pathways (MSigDB or your GMT)
   annotation gtf / ensembl    gene symbols and descriptions
                 │
           TSV tables
```

| Step | Command | R function | Output |
|---|---|---|---|
| Import | `data import` | `read_rsem_study()`, `write_study()` | study folder: counts, TPM, effective length, metadata |
| Summaries | `data summarize` | `summarize_study()` | per-sample and per-gene summaries |
| Sample structure | `expression pca`, `expression correlation` | `compute_pca()`, `compute_correlation()` | PCA scores, loadings, variance; sample × sample matrix |
| Normalized values | `expression normalize` | `normalize_expression()` | gene × sample matrix (`vst`, `normalized-counts`, `log2-tpm`) |
| Differential expression | `differential deseq2` | `differential_deseq2()` | one table per comparison, with the design matrix and contrasts |
| Pathways | `enrichment gsea`, `enrichment gsva` | `run_gsea()`, `run_gsva()` | GSEA results; pathway × sample scores |
| Annotation | `annotation gtf`, `annotation ensembl` | `annotate_gtf()`, `annotate_ensembl()` | one row per RSEM gene with symbol and status |

</details>

## Your data

<details>
<summary>Show details</summary>

```text
project/
├── data/
│   ├── S01.genes.results
│   └── S02.genes.results
└── metadata.tsv
```

```text
sample_id	condition	batch	age
S01	Control	A	51
S02	Drug	B	43
```

- **RSEM files:** the sample ID is the file name without `.genes.results`.
  Every file needs the standard RSEM columns and the same genes in the same
  order.
- **Metadata:** TSV or CSV with a `sample_id` column; every other column is
  optional. Every RSEM sample must be in the metadata and the other way round.
  `sample_id` is read as text (`001` stays `001`). A column of numbers becomes
  numeric and any other column categorical; override with
  `--p-categorical COLUMN` (for example time points `0, 6, 24` as groups) or
  `--p-numeric COLUMN`.

</details>

## A typical run

<details>
<summary>Show details</summary>

```bash
rsemflow data import --i-rsem-dir data/ --m-metadata metadata.tsv --o-study study/
rsemflow data summarize --i-study study/ --o-summary results/summary/
rsemflow expression pca --i-study study/ --o-pca results/pca/

rsemflow differential deseq2 --i-study study/ \
  --p-formula '~ batch + condition' --p-term condition \
  --p-reference condition::Control --o-differential results/differential/

rsemflow enrichment gsea \
  --i-differential results/differential/condition_Drug_vs_Control.tsv \
  --p-species human --p-collection hallmark --o-enrichment results/gsea.tsv
rsemflow enrichment gsva --i-study study/ --p-species human \
  --p-collection hallmark --o-scores results/gsva.tsv

rsemflow annotation gtf --i-study study/ --i-gtf genes.gtf.gz --o-annotation results/annotation.tsv
```

Option prefixes say what each option is: `--i-*` inputs, `--m-*` metadata,
`--p-*` parameters, `--o-*` outputs. Every level has `--help`
(`rsemflow differential deseq2 --help`). The
[Tutorial](https://daikikumakura.github.io/rsemflow/articles/tutorial.html)
runs every command on the example and shows what each one writes.

</details>

## Designs

<details>
<summary>Show details</summary>

`--p-formula` is the full model, `--p-term` chooses which results to write,
and `--p-reference COLUMN::LEVEL` sets reference levels.

| Design | Key options |
|---|---|
| Two groups | `--p-formula '~ condition' --p-term condition --p-reference condition::Control` |
| Batch-adjusted | `--p-formula '~ batch + condition' --p-term condition` |
| Paired | `--p-formula '~ patient + condition' --p-term condition` |
| Several groups | `--p-term treatment --p-reference treatment::Vehicle`; other pairs with `--p-contrast treatment::DrugB::DrugA` |
| Numeric covariate | `--p-formula '~ sex + age + score' --p-term score` |
| Interaction | `--p-formula '~ genotype * treatment' --p-term genotype:treatment` |
| Effect within one group | interaction formula + `--p-contrast treatment::Drug::Vehicle --p-at genotype::KO` |
| Any effect of time (LRT) | `--p-test lrt --p-formula '~ patient + time' --p-reduced '~ patient' --p-term time` |

Log fold changes are the ordinary DESeq2 estimates (no shrinkage). The design
matrix must be full rank; rsemflow does not judge whether you have enough
replicates. Every run also writes `design-matrix.tsv`, `contrasts.tsv` and
`analysis-info.tsv`. Each design is explained with its output in
[Designs](https://daikikumakura.github.io/rsemflow/articles/designs.html).

</details>

## Details that affect results

<details>
<summary>Show details</summary>

- **Gene sets:** `--p-collection` takes `hallmark`, `reactome`, `go-bp`,
  `go-mf`, `go-cc`, or a raw MSigDB code with `--p-subcollection`; human uses
  the human MSigDB and mouse the mouse MSigDB (`MH`, `M2`, `M5`). Use
  `--i-genesets file.gmt` for your own sets.
- **Gene matching:** by stable Ensembl ID (`ENSG00000141510.18` →
  `ENSG00000141510`). When several rows share an ID, GSEA keeps the largest
  absolute statistic and GSVA averages the expression.
- **Reproducibility:** GSEA estimates p-values by sampling with a fixed seed
  (`--p-seed`, default 42), so repeated runs give identical tables.
- **Short genes:** RSEM sets `effective_length` to 0 for genes shorter than
  the fragment length; rsemflow keeps the values and replaces those zeros with
  1 only inside DESeq2, as the DESeq2 vignette recommends.
- **PCA and correlation** use VST and the 500 most variable genes by default;
  GSVA uses `log2(TPM + 1)`.

</details>

## Validation

Besides the unit tests, rsemflow is checked on real RSEM output from ENCODE
(K562 cells, PTBP1 and SRSF1 shRNA knockdowns with a shared control; 59,526
genes). Each knockdown target is among the most significant decreases in its
own knockdown and nearly unchanged in the other. The script downloads the
data, runs the workflow and checks the results:
[`validation/encode-k562-knockdown/`](https://github.com/DaikiKumakura/rsemflow/tree/main/validation/encode-k562-knockdown).
R CMD check passes on Linux, macOS and Windows in CI.

## Scope

Included: RSEM gene-level import, expression summaries, annotation,
normalization, PCA, correlation, DESeq2, GSEA, GSVA/ssGSEA.

Not included: plots, reports, automatic sample removal or design selection,
log fold change shrinkage, alignment-level QC, transcript-level or splicing
analysis, tests on GSVA scores between groups.

## Related tools

<details>
<summary>Show details</summary>

- [rsemverse](https://daikikumakura.github.io/rsemverse/) reads the same RSEM
  files into a gene symbol by sample table in R and filters, normalizes and
  reshapes it for plotting.
- [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/) estimates
  cell-type composition from a TPM table.

</details>

## Development

<details>
<summary>Show details</summary>

```bash
make deps      # install dependencies
make document  # regenerate man/ and NAMESPACE
make test      # testthat suite
make check     # R CMD build + R CMD check
```

</details>

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210477>.
rsemflow coordinates other methods; please cite those you use as well
(`citation("rsemflow")` lists RSEM, tximport, DESeq2, fgsea, GSVA and MSigDB).

## License

MIT © Daiki Kumakura
