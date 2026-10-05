test_that("DESeq2 two-group differential analysis runs", {
  ex <- system.file("extdata", "example", package = "rsemflow")
  study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))

  x <- differential_deseq2(
    study,
    formula = ~ batch + condition,
    term = "condition",
    reference = "condition::Control"
  )
  expect_true("condition_Drug_vs_Control" %in% names(x$results))
  tab <- x$results[["condition_Drug_vs_Control"]]
  expect_true(all(c("gene_id", "log2FoldChange", "stat", "pvalue", "padj") %in% names(tab)))
})
