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

.check_geneset_overlap <- function(sets, ids, what) {
  universe <- unique(unlist(sets, use.names = FALSE))
  n_hit <- sum(ids %in% universe)
  if (!n_hit) {
    .stopf(
      paste0(
        "None of the %d %s match a gene in the gene sets. Gene sets use stable ",
        "Ensembl gene IDs (for example ENSG00000141510); check the species and ID type."
      ),
      length(ids), what
    )
  }
  n_hit
}

# Runs `expr` with a fixed RNG seed and restores the caller's RNG state.
.with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

#' Run preranked GSEA
#'
#' @param differential Differential-expression table or path.
#' @param species `auto`, `human`, or `mouse`.
#' @param collection MSigDB collection alias or raw collection code.
#' @param subcollection Optional MSigDB subcollection.
#' @param genesets Optional named gene-set list or GMT path.
#' @param rank Column used as ranking statistic, default `stat`.
#' @param min_size Minimum mapped gene-set size.
#' @param max_size Maximum mapped gene-set size.
#' @param eps Passed to `fgseaMultilevel`.
#' @param seed Random seed for `fgseaMultilevel`, which estimates p-values by
#'   sampling. The same seed gives the same table.
#' @return A data frame of GSEA results, sorted by `padj` and then by
#'   decreasing `|NES|`.
#' @details Gene IDs are matched to gene sets by stable Ensembl gene ID. When
#'   several rows share a stable ID, the row with the largest absolute ranking
#'   statistic is kept.
#'
#'   Unsigned statistics, such as the `stat` column of an LRT table, are
#'   rejected because preranked GSEA needs a direction.
#' @examples
#' ids <- sprintf("ENSG%011d", 1:100)
#' de <- data.frame(ensembl_gene_id = ids, stat = seq(5, -5, length.out = 100))
#' sets <- list(UP = ids[1:20], DOWN = ids[81:100], MIXED = ids[seq(1, 100, 5)])
#' run_gsea(de, genesets = sets, min_size = 10, max_size = 50)
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
    eps = 0,
    seed = 42L) {

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

  .check_geneset_overlap(gs$sets, names(stats_vec), "ranked genes")

  res <- .with_seed(as.integer(seed), fgsea::fgseaMultilevel(
    pathways = gs$sets,
    stats = stats_vec,
    minSize = as.integer(min_size),
    maxSize = as.integer(max_size),
    eps = as.numeric(eps),
    # Serial and without a progress bar: fgsea's own nproc = 1 switches to a
    # SnowParam worker with a progress bar on Windows.
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
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
#' Scores each gene set in each sample from `log2(TPM + 1)`. Rows sharing a
#' stable Ensembl ID are averaged first. GSVA uses a Gaussian kernel
#' (`kcdf = "Gaussian"`). Scores are descriptive; no test between groups is run.
#'
#' @param study A `rsemflow_study` or study path.
#' @param species `auto`, `human`, or `mouse`.
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
  .check_geneset_overlap(gs$sets, rownames(expr), "expression rows")

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
