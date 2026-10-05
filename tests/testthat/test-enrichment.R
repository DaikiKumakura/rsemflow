test_that("GSEA accepts custom Ensembl gene sets", {
  ids <- sprintf("ENSG%011d", 1:100)
  de <- data.frame(
    gene_id = paste0(ids, ".1"),
    ensembl_gene_id = ids,
    stat = seq(5, -5, length.out = 100),
    stringsAsFactors = FALSE
  )
  gs <- list(
    TOP = ids[1:20],
    BOTTOM = ids[81:100]
  )
  res <- run_gsea(de, genesets = gs, min_size = 10, max_size = 50)
  expect_true(all(c("pathway", "NES", "pval", "padj") %in% names(res)))
})

test_that("GSVA accepts custom Ensembl gene sets", {
  ex <- system.file("extdata", "example", package = "rsemflow")
  study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
  ids <- rsemflow:::.strip_ensembl_version(rownames(study$txi$abundance))
  gs <- list(
    SET_A = ids[1:25],
    SET_B = ids[26:50]
  )
  res <- run_gsva(study, genesets = gs, method = "gsva", min_size = 10, max_size = 50)
  expect_equal(nrow(res), 2)
  expect_true(all(study$sample_ids %in% names(res)))
})
