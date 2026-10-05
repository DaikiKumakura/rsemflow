test_that("metadata typing is predictable", {
  md <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "B", "A"),
    age = c(20, 30, 40),
    stringsAsFactors = FALSE
  )
  out <- rsemflow:::.coerce_metadata(md)
  expect_true(is.factor(out$group))
  expect_true(is.numeric(out$age))
})

test_that("rank deficient designs are rejected", {
  md <- data.frame(
    sample_id = paste0("S", 1:4),
    condition = factor(c("A", "A", "B", "B")),
    batch = factor(c("X", "X", "Y", "Y"))
  )
  expect_error(
    rsemflow:::.validate_design_matrix(md, ~ batch + condition),
    "not full rank"
  )
})
