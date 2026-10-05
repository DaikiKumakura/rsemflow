test_that("example RSEM data import works", {
  ex <- system.file("extdata", "example", package = "rsemflow")
  expect_true(nzchar(ex))
  study <- read_rsem_study(
    file.path(ex, "data"),
    file.path(ex, "metadata.tsv")
  )
  expect_s3_class(study, "rsemflow_study")
  expect_equal(length(study$sample_ids), 6)
  expect_equal(nrow(study$txi$counts), 200)
})
