.deduplicate_rank <- function(ids, stats) {
  keep <- !is.na(ids) & nzchar(ids) & !is.na(stats) & is.finite(stats)
  ids <- ids[keep]
  stats <- as.numeric(stats[keep])
  if (!length(stats)) .stopf("No finite ranked statistics are available for GSEA.")
  ord <- order(abs(stats), decreasing = TRUE)
  ids <- ids[ord]
  stats <- stats[ord]
  first <- !duplicated(ids)
  stats <- stats[first]
  names(stats) <- ids[first]
  sort(stats, decreasing = TRUE)
}

#' Run preranked GSEA
#'
#' @param differential Differential-expression table or path.
#' @param species `auto`, `human`, `mouse`, or `rat`.
#' @param collection MSigDB collection alias or raw collection code.
#' @param subcollection Optional MSigDB subcollection.
#' @param genesets Optional named gene-set list or GMT path.
#' @param rank Column used as ranking statistic, default `stat`.
#' @param min_size Minimum mapped gene-set size.
#' @param max_size Maximum mapped gene-set size.
#' @param eps Passed to `fgseaMultilevel`.
#' @return A data frame of GSEA results.
#' @export
run_gsea <- function(
    differential,
    species = "auto",
    collection = "hallmark",
    subcollection = NULL,
    genesets = NULL,
    rank = "stat",
    min_size = 15L,
    max_size = 500L,
    eps = 0) {

  de <- if (is.character(differential) && length(differential) == 1L) {
    .read_table_auto(differential)
  } else {
    as.data.frame(differential, stringsAsFactors = FALSE)
  }
  if (!rank %in% names(de)) .stopf("Ranking column '%s' is not present in the differential table.", rank)
  if ("test" %in% names(de) && all(toupper(as.character(de$test)) == "LRT") && identical(rank, "stat")) {
    .stopf(
      "LRT statistics are omnibus and unsigned, so they are not suitable for directional preranked GSEA. Use a Wald contrast/coefficient table or provide another signed ranking column."
    )
  }

  ids <- if ("ensembl_gene_id" %in% names(de)) {
    as.character(de$ensembl_gene_id)
  } else if ("gene_id" %in% names(de)) {
    .strip_ensembl_version(de$gene_id)
  } else {
    .stopf("Differential table must contain 'ensembl_gene_id' or 'gene_id'.")
  }
  stats_vec <- .deduplicate_rank(ids, de[[rank]])
  gs <- .resolve_genesets(
    genesets = genesets,
    species = species,
    collection = collection,
    subcollection = subcollection,
    ids = names(stats_vec)
  )

  res <- fgsea::fgseaMultilevel(
    pathways = gs$sets,
    stats = stats_vec,
    minSize = as.integer(min_size),
    maxSize = as.integer(max_size),
    eps = as.numeric(eps)
  )
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  if ("leadingEdge" %in% names(res)) {
    res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", character(1))
  }
  res$collection <- gs$collection
  res$subcollection <- gs$subcollection
  res$db_version <- paste(unique(gs$db_version), collapse = ",")
  if ("padj" %in% names(res)) {
    res <- res[order(res$padj, -abs(res$NES)), , drop = FALSE]
  }
  rownames(res) <- NULL
  res
}

#' Run GSVA or ssGSEA pathway scoring
#'
#' @param study A `rsemflow_study` or study path.
#' @param species `auto`, `human`, `mouse`, or `rat`.
#' @param collection MSigDB collection alias or raw collection code.
#' @param subcollection Optional MSigDB subcollection.
#' @param genesets Optional named list or GMT file.
#' @param method `gsva` or `ssgsea`.
#' @param min_size Minimum gene-set size.
#' @param max_size Maximum gene-set size.
#' @return A pathway-by-sample data frame with gene-set metadata columns.
#' @export
run_gsva <- function(
    study,
    species = "auto",
    collection = "hallmark",
    subcollection = NULL,
    genesets = NULL,
    method = c("gsva", "ssgsea"),
    min_size = 10L,
    max_size = 500L) {

  if (is.character(study) && length(study) == 1L) study <- read_study(study)
  method <- match.arg(method)

  ids <- .strip_ensembl_version(rownames(study$txi$abundance))
  expr <- log2(study$txi$abundance + 1)
  expr <- .collapse_matrix_by_id(expr, ids)

  gs <- .resolve_genesets(
    genesets = genesets,
    species = species,
    collection = collection,
    subcollection = subcollection,
    ids = rownames(expr)
  )

  param <- if (method == "gsva") {
    GSVA::gsvaParam(
      exprData = expr,
      geneSets = gs$sets,
      minSize = as.integer(min_size),
      maxSize = as.integer(max_size),
      kcdf = "Gaussian"
    )
  } else {
    GSVA::ssgseaParam(
      exprData = expr,
      geneSets = gs$sets,
      minSize = as.integer(min_size),
      maxSize = as.integer(max_size)
    )
  }

  scores <- GSVA::gsva(param, verbose = FALSE)
  if (!is.matrix(scores)) {
    scores <- as.matrix(scores)
  }
  out <- data.frame(
    pathway = rownames(scores),
    collection = gs$collection,
    subcollection = gs$subcollection,
    db_version = paste(unique(gs$db_version), collapse = ","),
    as.data.frame(scores, check.names = FALSE, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}
