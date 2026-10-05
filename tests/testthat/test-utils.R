test_that("only Ensembl gene versions are stripped", {
  x <- c("ENSG00000141510.18", "ENSMUSG00000059552.7", "ENSG00000182378.14_PAR_Y", "Y_RNA.1", "ABC")
  expect_equal(
    rsemflow:::.strip_ensembl_version(x),
    c("ENSG00000141510", "ENSMUSG00000059552", "ENSG00000182378.14_PAR_Y", "Y_RNA.1", "ABC")
  )
})

test_that("species are inferred from Ensembl prefixes", {
  expect_equal(rsemflow:::.resolve_species("auto", "ENSG00000141510.1"), "human")
  expect_equal(rsemflow:::.resolve_species("auto", "ENSMUSG00000059552"), "mouse")
  expect_equal(rsemflow:::.resolve_species("Rattus norvegicus"), "rat")
  expect_error(rsemflow:::.resolve_species("auto", "geneA"), "could not be inferred")
  expect_error(rsemflow:::.resolve_species("zebrafish"), "Unsupported species")
})

test_that("collection aliases resolve", {
  x <- rsemflow:::.collection_alias("reactome")
  expect_equal(x$collection, "C2")
  expect_equal(x$subcollection, "CP:REACTOME")
  m <- rsemflow:::.collection_alias("hallmark", db_species = "MM")
  expect_equal(m$collection, "MH")
  expect_equal(rsemflow:::.collection_alias("go-bp", db_species = "MM")$collection, "M5")
  y <- rsemflow:::.collection_alias("C7", "IMMUNESIGDB")
  expect_equal(y, list(collection = "C7", subcollection = "IMMUNESIGDB"))
})

test_that("duplicated IDs are averaged", {
  m <- matrix(c(1, 3, 5, 2, 4, 6), 3, dimnames = list(NULL, c("a", "b")))
  out <- rsemflow:::.collapse_matrix_by_id(m, c("x", "x", "y"))
  expect_equal(out["x", ], c(a = 2, b = 3))
  expect_equal(out["y", ], c(a = 5, b = 6))
})

test_that("GMT files are parsed", {
  f <- tempfile(fileext = ".gmt")
  writeLines(c("SET1\tdesc\tG1\tG2\tG2", "EMPTY\tdesc", "SET2\thttp://x\tG3"), f)
  sets <- rsemflow:::.read_gmt(f)
  expect_equal(sets, list(SET1 = c("G1", "G2"), SET2 = "G3"))
})
