test_that("PCA and correlation return tables", {
  ex <- system.file("extdata", "example", package = "rsemflow")
  study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))

  p <- compute_pca(study, ntop = 50, components = 3)
  expect_equal(nrow(p$scores), 6)
  expect_true(all(c("PC1", "PC2") %in% names(p$scores)))
  expect_equal(nrow(p$loadings), 50)

  cr <- compute_correlation(study)
  expect_equal(dim(cr), c(6, 6))
  expect_equal(diag(cr), rep(1, 6), tolerance = 1e-8)
})
