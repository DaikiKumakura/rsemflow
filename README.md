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
head(de$results$condition_Drug_vs_Control)
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
