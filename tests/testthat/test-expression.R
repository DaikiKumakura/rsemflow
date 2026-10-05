test_that("normalization methods return gene-by-sample matrices", {
  study <- example_study()
  for (m in c("vst", "normalized-counts", "log2-tpm")) {
    x <- normalize_expression(study, m)
    expect_equal(dim(x), c(200L, 6L))
    expect_equal(colnames(x), study$sample_ids)
  }
  expect_equal(
    normalize_expression(study, "log2-tpm"),
    log2(study$txi$abundance + 1)
  )
})

test_that("summaries describe each sample and gene", {
  s <- summarize_study(example_study())
  expect_equal(nrow(s$samples), 6L)
  expect_equal(nrow(s$genes), 200L)
  expect_equal(s$samples$detected_genes + s$samples$undetected_genes, rep(200, 6))
})

test_that("PCA and correlation return tables", {
  study <- example_study()
  p <- compute_pca(study, ntop = 50, components = 3)
  expect_equal(nrow(p$scores), 6)
  expect_equal(names(p$scores), c("sample_id", "PC1", "PC2", "PC3"))
  expect_equal(nrow(p$loadings), 50)
  expect_equal(sum(p$variance$variance_fraction), 1)
  expect_equal(p$ntop, 50L)

  cr <- compute_correlation(study)
  expect_equal(dim(cr), c(6L, 6L))
  expect_equal(unname(diag(cr)), rep(1, 6), tolerance = 1e-8)
  sp <- compute_correlation(study, transform = "log2-tpm", method = "spearman")
  expect_true(isSymmetric(sp))
})
