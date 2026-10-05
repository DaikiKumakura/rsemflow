#' rsemflow: metadata-driven downstream analysis for RSEM
#'
#' rsemflow imports RSEM `*.genes.results` files and ordinary sample metadata,
#' then produces tabular downstream-analysis results: expression summaries,
#' normalization, PCA, sample correlation, DESeq2 differential expression,
#' preranked GSEA, GSVA/ssGSEA scores, and gene annotation.
#'
#' The same functions back the `rsemflow` command-line tool (see
#' [cli_main()]). The package does not draw figures or write reports; every
#' output is a tab-separated table for use in any downstream tool.
#'
#' @keywords internal
"_PACKAGE"
