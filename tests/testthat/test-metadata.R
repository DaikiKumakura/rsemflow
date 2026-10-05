test_that("metadata typing is predictable", {
  md <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "B", "A"),
    age = c(20, 30, 40),
    time = c("0", "6", "24"),
    stringsAsFactors = FALSE
  )
  out <- rsemflow:::.coerce_metadata(md)
  expect_true(is.factor(out$group))
  expect_true(is.numeric(out$age))
  expect_true(is.numeric(out$time))

  forced <- rsemflow:::.coerce_metadata(md, categorical = "time")
  expect_true(is.factor(forced$time))
  expect_error(rsemflow:::.coerce_metadata(md, numeric = "group"), "cannot be converted")
  expect_error(rsemflow:::.coerce_metadata(md, categorical = "nope"), "Unknown --p-categorical")
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

test_that("missing values in model variables are rejected", {
  md <- data.frame(sample_id = paste0("S", 1:4), x = c(1, NA, 3, 4))
  expect_error(rsemflow:::.validate_formula_metadata(md, ~x), "Missing metadata values")
})

test_that("reference levels are applied and validated", {
  md <- data.frame(sample_id = paste0("S", 1:4), g = factor(c("b", "a", "b", "a")))
  out <- rsemflow:::.apply_reference_levels(md, "g::b")
  expect_equal(levels(out$g)[1], "b")
  expect_error(rsemflow:::.apply_reference_levels(md, "g::z"), "not present")
  expect_error(rsemflow:::.apply_reference_levels(md, "g"), "expects 2 fields")
})

test_that("duplicated or empty sample IDs are rejected", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tg", "A\tx", "A\ty"), f)
  expect_error(rsemflow:::.read_metadata(f), "duplicated sample_id")
  writeLines(c("sample_id\tg", "A\tx", "\ty"), f)
  expect_error(rsemflow:::.read_metadata(f), "missing or empty")
})

test_that("CSV metadata is read", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("sample_id,g", "A,x", "B,y"), f)
  md <- rsemflow:::.read_metadata(f)
  expect_equal(md$g, c("x", "y"))
})
