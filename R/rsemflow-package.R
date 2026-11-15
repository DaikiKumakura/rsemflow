#' rsemflow: standard bulk RNA-seq results from RSEM files and a sample table
#'
#' rsemflow takes RSEM `*.genes.results` files and an ordinary sample
#' metadata table and produces the usual downstream results as tables. The
#' same functions back the `rsemflow` command-line tool ([cli_main()]).
#'
#' @section Workflow:
#' 1. Import: [read_rsem_study()] checks the files against the metadata;
#'    [write_study()] and [read_study()] store and reload the study.
#' 2. Summaries: [summarize_study()].
#' 3. Sample structure: [compute_pca()], [compute_correlation()].
#' 4. Normalized values: [normalize_expression()] (VST, normalized counts,
#'    log2 TPM).
#' 5. Differential expression: [differential_deseq2()], with the design
#'    stated explicitly.
#' 6. Pathways: [run_gsea()], [run_gsva()].
#' 7. Annotation: [annotate_gtf()], [annotate_ensembl()].
#'
#' The design is always given by the user, never guessed, and inconsistent
#' input stops the run with a message. The package writes tables only; it does
#' not draw figures or write reports.
#'
#' @keywords internal
"_PACKAGE"
