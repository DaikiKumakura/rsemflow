.make_dds_for_transform <- function(study) {
  md <- study$metadata
  DESeq2::DESeqDataSetFromTximport(
    study$txi,
    colData = md,
    design = stats::as.formula("~ 1")
  )
}

#' Normalize or transform expression values
#'
#' @param study A `rsemflow_study` or study path.
#' @param method One of `vst`, `normalized-counts`, or `log2-tpm`.
#' @return A gene-by-sample numeric matrix.
#' @export
normalize_expression <- function(study, method = c("vst", "normalized-counts", "log2-tpm")) {
  if (is.character(study)) study <- read_study(study)
  method <- match.arg(method)

  if (method == "log2-tpm") {
    out <- log2(study$txi$abundance + 1)
    return(out)
  }

  dds <- .make_dds_for_transform(study)
  dds <- DESeq2::estimateSizeFactors(dds)

  if (method == "normalized-counts") {
    return(DESeq2::counts(dds, normalized = TRUE))
  }

  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  SummarizedExperiment::assay(vsd)
}

#' Compute PCA tables from transformed expression
#'
#' @param study A `rsemflow_study` or path.
#' @param transform `vst` or `log2-tpm`.
#' @param ntop Number of most variable genes to use.
#' @param components Maximum number of PCs to return.
#' @return A list containing `scores`, `loadings`, and `variance`.
#' @export
compute_pca <- function(study, transform = c("vst", "log2-tpm"), ntop = 500L, components = 10L) {
  if (is.character(study)) study <- read_study(study)
  transform <- match.arg(transform)
  ntop <- as.integer(ntop)
  components <- as.integer(components)
  if (ntop < 2L) .stopf("ntop must be >= 2.")
  if (components < 1L) .stopf("components must be >= 1.")

  mat <- normalize_expression(study, method = transform)
  rv <- apply(mat, 1L, stats::var, na.rm = TRUE)
  rv[is.na(rv)] <- -Inf
  ord <- order(rv, decreasing = TRUE)
  keep <- ord[seq_len(min(ntop, length(ord)))]
  selected <- mat[keep, , drop = FALSE]

  p <- stats::prcomp(t(selected), center = TRUE, scale. = FALSE)
  k <- min(components, .available_pc_count(p))
  pcs <- seq_len(k)
  pc_names <- paste0("PC", pcs)

  scores <- data.frame(
    sample_id = rownames(p$x),
    p$x[, pcs, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(scores)[-1] <- pc_names

  loadings <- data.frame(
    gene_id = rownames(p$rotation),
    ensembl_gene_id = .strip_ensembl_version(rownames(p$rotation)),
    p$rotation[, pcs, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(loadings)[-(1:2)] <- pc_names

  eigen <- p$sdev^2
  frac <- eigen / sum(eigen)
  variance <- data.frame(
    component = paste0("PC", seq_along(eigen)),
    eigenvalue = eigen,
    variance_fraction = frac,
    cumulative_variance = cumsum(frac),
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    loadings = loadings,
    variance = variance,
    transform = transform,
    ntop = min(ntop, nrow(mat))
  )
}

#' Compute sample correlation
#'
#' @param study A `rsemflow_study` or study path.
#' @param transform `vst` or `log2-tpm`.
#' @param method `pearson` or `spearman`.
#' @return A sample-by-sample correlation matrix.
#' @export
compute_correlation <- function(study, transform = c("vst", "log2-tpm"), method = c("pearson", "spearman")) {
  if (is.character(study)) study <- read_study(study)
  transform <- match.arg(transform)
  method <- match.arg(method)
  mat <- normalize_expression(study, method = transform)
  stats::cor(mat, method = method, use = "pairwise.complete.obs")
}
