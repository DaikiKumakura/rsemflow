test_that("Ensembl version suffixes are stripped", {
  x <- c("ENSG00000141510.18", "ENSMUSG00000059552.7", "ABC")
  expect_equal(
    rsemflow:::.strip_ensembl_version(x),
    c("ENSG00000141510", "ENSMUSG00000059552", "ABC")
  )
})

test_that("collection aliases resolve", {
  x <- rsemflow:::.collection_alias("reactome")
  expect_equal(x$collection, "C2")
  expect_equal(x$subcollection, "CP:REACTOME")
})
