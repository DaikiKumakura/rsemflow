test_that("example RSEM data import works", {
  study <- example_study()
  expect_s3_class(study, "rsemflow_study")
  expect_equal(study$sample_ids, c("C1", "C2", "C3", "D1", "D2", "D3"))
  expect_equal(dim(study$txi$counts), c(200L, 6L))
  expect_true(is.factor(study$metadata$condition))
  expect_true(is.numeric(study$metadata$age))
  expect_output(print(study), "samples:  6")
})

test_that("write_study and read_study round-trip", {
  study <- example_study()
  out <- test_tempdir()
  write_study(study, out)
  expect_setequal(
    list.files(out),
    c("study.rds", "metadata.tsv", "samples.tsv", "expected_count.tsv",
      "tpm.tsv", "effective_length.tsv", "study-info.tsv")
  )
  back <- read_study(out)
  expect_identical(back$txi$counts, study$txi$counts)
  counts_tsv <- utils::read.delim(file.path(out, "expected_count.tsv"), check.names = FALSE)
  expect_equal(names(counts_tsv), c("gene_id", study$sample_ids))
  expect_equal(nrow(counts_tsv), 200L)
})

test_that("read_study rejects non-study objects", {
  path <- tempfile(fileext = ".rds")
  saveRDS(list(a = 1), path)
  expect_error(read_study(path), "not an RSEMflow study")
})

test_that("sample mismatches between files and metadata stop the import", {
  counts <- simulate_counts(20, c("A", "B", "C"))
  md <- data.frame(sample_id = c("A", "B", "X"), group = c("g1", "g1", "g2"))
  dir <- write_rsem_dataset(counts, md)
  expect_error(
    read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv")),
    "RSEM samples missing from metadata: C"
  )
})

test_that("gene tables must match across samples", {
  counts <- simulate_counts(20, c("A", "B"))
  md <- data.frame(sample_id = c("A", "B"), group = c("g1", "g2"))
  dir <- write_rsem_dataset(counts, md)
  f <- file.path(dir, "data", "B.genes.results")
  tab <- utils::read.delim(f, check.names = FALSE)
  utils::write.table(tab[rev(seq_len(nrow(tab))), ], f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(
    read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv")),
    "differ between RSEM files"
  )
})

test_that("missing RSEM columns are reported", {
  counts <- simulate_counts(5, c("A", "B"))
  md <- data.frame(sample_id = c("A", "B"), group = c("g1", "g2"))
  dir <- write_rsem_dataset(counts, md)
  f <- file.path(dir, "data", "A.genes.results")
  tab <- utils::read.delim(f, check.names = FALSE)
  tab$FPKM <- NULL
  utils::write.table(tab, f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(
    read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv")),
    "missing required columns: FPKM"
  )
})

test_that("numeric-looking sample IDs keep leading zeros", {
  counts <- simulate_counts(20, c("001", "002", "003", "004"))
  md <- data.frame(sample_id = c("001", "002", "003", "004"), group = c("a", "a", "b", "b"))
  dir <- write_rsem_dataset(counts, md)
  study <- read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv"))
  expect_equal(study$metadata$sample_id, c("001", "002", "003", "004"))
})

test_that("genes with zero effective length do not break DESeq2 steps", {
  ids <- paste0("S", 1:6)
  counts <- simulate_counts(60, ids)
  counts[1:3, ] <- 0
  len <- c(0, 0, 0, rep(1000, 57))
  md <- data.frame(sample_id = ids, group = rep(c("a", "b"), each = 3))
  dir <- write_rsem_dataset(counts, md, effective_length = len)
  study <- read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv"))
  expect_true(any(study$txi$length == 0))

  vst <- normalize_expression(study, "vst")
  expect_equal(dim(vst), c(60L, 6L))
  de <- differential_deseq2(study, ~group, "group")
  expect_equal(nrow(de$results$group_b_vs_a), 60L)

  out <- test_tempdir()
  write_study(study, out)
  info <- utils::read.delim(file.path(out, "study-info.tsv"))
  expect_equal(info$value[info$key == "zero_effective_length_values"], "18")
})
