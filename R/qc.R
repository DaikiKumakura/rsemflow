#' Summarize expression-level QC statistics
#'
#' @param study A `rsemflow_study` or path to a study.
#' @return A list with `samples` and `genes` data frames.
#' @export
summarize_study <- function(study) {
  if (is.character(study)) study <- read_study(study)
  counts <- study$txi$counts
  tpm <- study$txi$abundance

  sample_tbl <- data.frame(
    sample_id = colnames(counts),
    total_expected_count = colSums(counts, na.rm = TRUE),
    detected_genes = colSums(counts > 0, na.rm = TRUE),
    undetected_genes = colSums(counts <= 0 | is.na(counts), na.rm = TRUE),
    detected_fraction = colMeans(counts > 0, na.rm = TRUE),
    tpm_q25 = apply(tpm, 2L, stats::quantile, probs = 0.25, na.rm = TRUE, names = FALSE),
    tpm_median = apply(tpm, 2L, stats::median, na.rm = TRUE),
    tpm_q75 = apply(tpm, 2L, stats::quantile, probs = 0.75, na.rm = TRUE, names = FALSE),
    tpm_sum = colSums(tpm, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  gene_tbl <- data.frame(
    gene_id = rownames(counts),
    ensembl_gene_id = .strip_ensembl_version(rownames(counts)),
    mean_count = rowMeans(counts, na.rm = TRUE),
    median_count = apply(counts, 1L, stats::median, na.rm = TRUE),
    mean_tpm = rowMeans(tpm, na.rm = TRUE),
    median_tpm = apply(tpm, 1L, stats::median, na.rm = TRUE),
    detected_samples = rowSums(counts > 0, na.rm = TRUE),
    detected_fraction = rowMeans(counts > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(samples = sample_tbl, genes = gene_tbl)
}
