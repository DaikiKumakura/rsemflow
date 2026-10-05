ranked_table <- function() {
  ids <- sprintf("ENSG%011d", 1:100)
  data.frame(
    gene_id = paste0(ids, ".1"),
    ensembl_gene_id = ids,
    stat = seq(5, -5, length.out = 100),
    test = "Wald",
    stringsAsFactors = FALSE
  )
}

test_that("GSEA accepts custom gene sets and finds the simulated direction", {
  de <- ranked_table()
  gs <- list(TOP = de$ensembl_gene_id[1:20], BOTTOM = de$ensembl_gene_id[81:100])
  res <- run_gsea(de, genesets = gs, min_size = 10, max_size = 50)
  expect_true(all(c("pathway", "NES", "pval", "padj", "leadingEdge", "collection") %in% names(res)))
  expect_gt(res$NES[res$pathway == "TOP"], 0)
  expect_lt(res$NES[res$pathway == "BOTTOM"], 0)
  expect_type(res$leadingEdge, "character")
  expect_equal(unique(res$collection), "custom")
})

test_that("GSEA is reproducible with a fixed seed", {
  de <- ranked_table()
  gs <- list(A = de$ensembl_gene_id[seq(1, 100, 3)], B = de$ensembl_gene_id[seq(2, 100, 4)])
  expect_identical(
    run_gsea(de, genesets = gs, min_size = 5, seed = 7),
    run_gsea(de, genesets = gs, min_size = 5, seed = 7)
  )
})

test_that("GSEA leaves the caller's random stream untouched", {
  de <- ranked_table()
  gs <- list(A = de$ensembl_gene_id[seq(1, 100, 3)])
  set.seed(99)
  expected <- stats::runif(1)
  set.seed(99)
  run_gsea(de, genesets = gs, min_size = 5)
  expect_equal(stats::runif(1), expected)
})

test_that("duplicated stable IDs keep the largest absolute statistic", {
  stats <- rsemflow:::.deduplicate_rank(c("g1", "g1", "g2"), c(1, -3, 2))
  expect_equal(stats, c(g2 = 2, g1 = -3))
})

test_that("LRT tables and unmatched IDs are rejected", {
  de <- ranked_table()
  gs <- list(TOP = de$ensembl_gene_id[1:20])
  lrt <- de
  lrt$test <- "LRT"
  expect_error(run_gsea(lrt, genesets = gs), "unsigned")
  expect_error(run_gsea(de, genesets = list(X = paste0("GENE", 1:20))), "None of the")
  expect_error(run_gsea(de, genesets = gs, rank = "nope"), "not present")
})

test_that("GSVA and ssGSEA accept custom gene sets", {
  study <- example_study()
  ids <- rsemflow:::.strip_ensembl_version(rownames(study$txi$abundance))
  gs <- list(SET_A = ids[1:25], SET_B = ids[26:50])
  for (m in c("gsva", "ssgsea")) {
    res <- run_gsva(study, genesets = gs, method = m, min_size = 10, max_size = 50)
    expect_equal(res$pathway, c("SET_A", "SET_B"))
    expect_equal(
      names(res),
      c("pathway", "collection", "subcollection", "db_version", study$sample_ids)
    )
  }
})

skip_if_no_msigdb <- function() {
  ok <- tryCatch(nrow(msigdbr::msigdbr_collections()) > 0L, error = function(e) FALSE)
  skip_if_not(ok, "msigdbr data unavailable")
}

test_that("MSigDB hallmark sets load through msigdbr", {
  skip_on_cran()
  skip_if_no_msigdb()
  gs <- rsemflow:::.get_msigdbr_genesets("human", "hallmark")
  expect_equal(length(gs$sets), 50L)
  expect_true(all(grepl("^HALLMARK_", names(gs$sets))))
  expect_true(all(grepl("^ENSG", unlist(gs$sets[1:3]))))
})

test_that("mouse aliases resolve to the mouse MSigDB collections", {
  skip_on_cran()
  skip_if_no_msigdb()
  gs <- rsemflow:::.get_msigdbr_genesets("mouse", "hallmark")
  expect_equal(gs$collection, "MH")
  expect_true(all(grepl("^ENSMUSG", unlist(gs$sets[1:3]))))
})
