test_that("help is printed at every level", {
  expect_output(cli_main(character()), "Modules and actions")
  expect_output(cli_main("--version"), as.character(utils::packageVersion("rsemflow")))
  expect_output(cli_main("differential"), "deseq2")
  expect_output(cli_main(c("differential", "deseq2", "--help")), "--p-formula")
  expect_output(cli_main(c("enrichment", "gsea", "-h")), "--p-seed")
})

test_that("bad command lines are rejected with clear messages", {
  expect_error(cli_main(c("data", "nope")), "Unknown module/action")
  expect_error(
    cli_main(c("data", "import", "--i-rsem-dir", "x")),
    "Missing required option --m-metadata"
  )
  expect_error(cli_main(c("data", "import", "--bogus", "x")), "Unknown option")
  expect_error(cli_main(c("data", "summarize", "--i-study")), "requires a value")
  expect_error(
    cli_main(c("data", "summarize", "--i-study", "a", "--i-study", "b", "--o-summary", "o")),
    "only once"
  )
  expect_error(
    cli_main(c("expression", "pca", "--i-study", "s", "--o-pca", "o", "--p-ntop", "many")),
    "expects an integer"
  )
})

test_that("the CLI runs the example workflow end to end", {
  ex <- example_dir()
  out <- test_tempdir()
  study <- file.path(out, "study")
  run <- function(...) cli_main(c(...))

  expect_output(
    run("data", "import",
        "--i-rsem-dir", file.path(ex, "data"),
        "--m-metadata", file.path(ex, "metadata.tsv"),
        "--o-study", study),
    "Imported 6 samples and 200 genes"
  )
  expect_output(run("data", "summarize", "--i-study", study, "--o-summary", file.path(out, "summary")))
  expect_output(run("expression", "normalize", "--i-study", study, "--p-method", "log2-tpm",
                    "--o-expression", file.path(out, "log2tpm.tsv")))
  expect_output(run("expression", "pca", "--i-study", study, "--p-ntop", "100",
                    "--o-pca", file.path(out, "pca")))
  expect_output(run("expression", "correlation", "--i-study", study,
                    "--o-correlation", file.path(out, "correlation.tsv")))
  expect_output(
    run("differential", "deseq2", "--i-study", study,
        "--p-formula", "~ batch + condition", "--p-term", "condition",
        "--p-reference", "condition::Control",
        "--o-differential", file.path(out, "differential")),
    "Wrote 1 differential result table"
  )

  expect_setequal(list.files(file.path(out, "summary")), c("samples.tsv", "genes.tsv"))
  expect_setequal(
    list.files(file.path(out, "pca")),
    c("scores.tsv", "loadings.tsv", "variance.tsv", "analysis-info.tsv")
  )
  expect_setequal(
    list.files(file.path(out, "differential")),
    c("design-matrix.tsv", "contrasts.tsv", "analysis-info.tsv", "condition_Drug_vs_Control.tsv")
  )
  info <- utils::read.delim(file.path(out, "differential", "analysis-info.tsv"))
  expect_equal(info$value[info$key == "formula"], "~batch + condition")
  expect_equal(info$value[info$key == "alpha"], "0.05")
  de_path <- file.path(out, "differential", "condition_Drug_vs_Control.tsv")
  de <- utils::read.delim(de_path)
  expect_equal(nrow(de), 200L)

  ids <- unique(de$ensembl_gene_id)
  gmt <- file.path(out, "sets.gmt")
  writeLines(c(
    paste(c("SET_A", "na", ids[1:30]), collapse = "\t"),
    paste(c("SET_B", "na", ids[31:60]), collapse = "\t")
  ), gmt)
  expect_output(
    run("enrichment", "gsea", "--i-differential", de_path, "--i-genesets", gmt,
        "--p-min-size", "10", "--o-enrichment", file.path(out, "gsea.tsv")),
    "Wrote 2 GSEA results"
  )
  expect_output(
    run("enrichment", "gsva", "--i-study", study, "--i-genesets", gmt,
        "--p-method", "ssgsea", "--o-scores", file.path(out, "gsva.tsv")),
    "Wrote 2 pathway score rows"
  )
})
