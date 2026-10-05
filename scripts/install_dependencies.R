cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(cran_repo) || identical(cran_repo, "@CRAN@")) {
  cran_repo <- "https://cloud.r-project.org"
}

cran <- c("data.table", "msigdbr", "testthat")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = cran_repo)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = cran_repo)
}

bioc <- c(
  "tximport",
  "DESeq2",
  "SummarizedExperiment",
  "AnnotationDbi",
  "fgsea",
  "GSVA",
  "org.Hs.eg.db",
  "org.Mm.eg.db"
)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
