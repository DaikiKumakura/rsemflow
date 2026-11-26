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
res <- de$results$condition_Drug_vs_Control
head(res[, c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")], 4)
```

```text
            gene_id baseMean log2FoldChange       pvalue         padj
1 ENSG00000000001.2 169.1780       2.397434 3.199762e-05 0.0007999406
2 ENSG00000000002.3 208.8794       2.232664 2.306482e-04 0.0038441366
3 ENSG00000000003.1 248.5807       2.127022 1.475414e-03 0.0126338835
4 ENSG00000000004.2 288.2820       2.053341 7.830909e-03 0.0522060612
```

Each row is a gene: its mean expression, the log2 fold change of Drug against Control (adjusted for `batch`),
and the p-value and the adjusted p-value of the DESeq2 Wald test. `summary` holds the per-sample and per-gene
summaries and `pca` the principal components of the study.

The same from the command line, writing every table to a folder:

```bash
bash inst/extdata/example/run_example.sh example-output
```

## Documentation

[Package website](https://daikikumakura.github.io/rsemflow/) · [Function reference](https://daikikumakura.github.io/rsemflow/reference/) · [Changelog](https://daikikumakura.github.io/rsemflow/news/) · [Releases](https://github.com/DaikiKumakura/rsemflow/releases) · [R-universe](https://daikikumakura.r-universe.dev/rsemflow) · [Usage guide](https://github.com/DaikiKumakura/rsemflow/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Ran the analysis on ENCODE K562 RSEM output for PTBP1 and SRSF1 knockdowns and ran package checks in CI.
- **Results:** Each target was among the most downregulated genes in its own knockdown and nearly unchanged in the other.
- **Limitations:** The demonstration covers one experimental design and does not validate every design or gene-set analysis.
- **Reproduction:** [validation/encode-k562-knockdown](https://github.com/DaikiKumakura/rsemflow/tree/main/validation/encode-k562-knockdown).

## Citation

Kumakura, D. (2026). *rsemflow: Metadata-Driven Downstream Analysis for RSEM Gene-Level Results* [Software]. [DOI: 10.5281/zenodo.24210477](https://doi.org/10.5281/zenodo.24210477).

Citation metadata: [CITATION.cff](CITATION.cff).

Also cite the methods used. Run `citation("rsemflow")` for the list.

## License

This project is licensed under MIT. See [LICENSE.md](LICENSE.md) for the license terms.
