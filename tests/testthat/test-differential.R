two_group_dataset <- function() {
  ids <- c(paste0("C", 1:3), paste0("D", 1:3))
  effect <- matrix(0, 200, 6)
  effect[1:20, 4:6] <- 2
  counts <- simulate_counts(200, ids, effect, seed = 3L)
  md <- data.frame(sample_id = ids, condition = rep(c("Control", "Drug"), each = 3))
  dir <- write_rsem_dataset(counts, md)
  read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv"))
}

test_that("two-group analysis recovers the simulated direction", {
  study <- two_group_dataset()
  x <- differential_deseq2(study, ~condition, "condition", reference = "condition::Control")
  expect_named(x$results, "condition_Drug_vs_Control")
  tab <- x$results$condition_Drug_vs_Control
  expect_equal(
    names(tab)[1:6],
    c("gene_id", "ensembl_gene_id", "term", "contrast", "coefficient", "test")
  )
  expect_true(all(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in% names(tab)))
  expect_gt(stats::median(tab$log2FoldChange[1:20]), 1.5)
  expect_lt(abs(stats::median(tab$log2FoldChange[21:200])), 0.3)
  expect_equal(x$manifest$contrast, "condition::Drug::Control")
})

test_that("the example batch-adjusted analysis runs", {
  x <- differential_deseq2(
    example_study(), ~ batch + condition, "condition",
    reference = "condition::Control"
  )
  expect_named(x$results, "condition_Drug_vs_Control")
  expect_equal(ncol(x$design_matrix), 3L)
})

test_that("multi-group terms give one table per level plus explicit contrasts", {
  ids <- sprintf("S%02d", 1:9)
  md <- data.frame(sample_id = ids, treatment = rep(c("Vehicle", "DrugA", "DrugB"), each = 3))
  dir <- write_rsem_dataset(simulate_counts(100, ids), md)
  study <- read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv"))
  x <- differential_deseq2(
    study, ~treatment, "treatment",
    reference = "treatment::Vehicle",
    contrast = c("treatment::DrugB::DrugA", "treatment::DrugA::Vehicle")
  )
  expect_setequal(
    names(x$results),
    c("treatment_DrugA_vs_Vehicle", "treatment_DrugB_vs_Vehicle", "treatment_DrugB_vs_DrugA")
  )
  expect_equal(nrow(x$manifest), 3L)
  expect_error(
    differential_deseq2(study, ~treatment, "treatment", contrast = "treatment::DrugC::Vehicle"),
    "Unknown level"
  )
})

test_that("numeric terms report the per-unit coefficient", {
  x <- differential_deseq2(example_study(), ~ age + condition, "age")
  expect_named(x$results, "age")
  expect_equal(x$manifest$contrast, "per_unit")
  expect_equal(x$manifest$coefficient, "age")
})

test_that("interaction terms and simple effects use the right coefficients", {
  dir <- interaction_dataset()
  study <- read_rsem_study(file.path(dir, "data"), file.path(dir, "metadata.tsv"))
  refs <- c("genotype::WT", "treatment::Vehicle")

  inter <- differential_deseq2(study, ~ genotype * treatment, "genotype:treatment", reference = refs)
  expect_equal(inter$manifest$coefficient, "genotypeKO.treatmentDrug")
  expect_gt(stats::median(inter$results[[1]]$log2FoldChange[1:50]), 1.5)

  at_ko <- differential_deseq2(
    study, ~ genotype * treatment, "treatment",
    reference = refs, contrast = "treatment::Drug::Vehicle", at = "genotype::KO"
  )
  expect_named(at_ko$results, "treatment_Drug_vs_Vehicle_at_genotype_KO")
  expect_equal(at_ko$manifest$coefficient, "treatment_Drug_vs_Vehicle+genotypeKO.treatmentDrug")

  # The same quantity from the reparameterized model with KO as reference.
  direct <- suppressWarnings(differential_deseq2(
    study, ~ genotype * treatment, "treatment",
    reference = c("genotype::KO", "treatment::Vehicle")
  ))
  expect_equal(
    at_ko$results[[1]]$log2FoldChange,
    direct$results$treatment_Drug_vs_Vehicle$log2FoldChange,
    tolerance = 1e-3
  )

  at_wt <- differential_deseq2(
    study, ~ genotype * treatment, "treatment",
    reference = refs, contrast = "treatment::Drug::Vehicle", at = "genotype::WT"
  )
  expect_lt(abs(stats::median(at_wt$results[[1]]$log2FoldChange[1:50])), 0.5)

  expect_warning(
    differential_deseq2(study, ~ genotype * treatment, "treatment", reference = refs),
    "interacts with genotype"
  )
  expect_error(
    differential_deseq2(
      study, ~ genotype + treatment, "treatment",
      reference = refs, contrast = "treatment::Drug::Vehicle", at = "genotype::KO"
    ),
    "not part of the formula"
  )
})

test_that("LRT returns an omnibus table without fold changes", {
  ids <- sprintf("S%02d", 1:9)
  md <- data.frame(
    sample_id = ids,
    patient = rep(c("P1", "P2", "P3"), 3),
    time = rep(c(0, 6, 24), each = 3)
  )
  dir <- write_rsem_dataset(simulate_counts(100, ids), md)
  study <- read_rsem_study(
    file.path(dir, "data"), file.path(dir, "metadata.tsv"),
    categorical = "time"
  )
  expect_true(is.factor(study$metadata$time))
  x <- differential_deseq2(
    study, ~ patient + time, "time",
    test = "lrt", reduced = ~patient
  )
  expect_named(x$results, "time_omnibus")
  tab <- x$results$time_omnibus
  expect_true(all(is.na(tab$log2FoldChange)))
  expect_true(all(tab$test == "LRT"))

  expect_error(
    differential_deseq2(study, ~ patient + time, "time", test = "lrt"),
    "--p-reduced is required"
  )
  expect_error(
    differential_deseq2(study, ~ patient + time, "time", test = "lrt", reduced = ~time),
    "must be removed"
  )
})

test_that("invalid requests fail before fitting", {
  study <- example_study()
  expect_error(differential_deseq2(study, ~condition, "batch"), "not used in the formula")
  expect_error(differential_deseq2(study, ~condition, "nope"), "unknown metadata columns")
  expect_error(differential_deseq2(study, ~condition, "condition", alpha = 2), "alpha")
})
