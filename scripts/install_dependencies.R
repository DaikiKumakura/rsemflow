# Install the CRAN and Bioconductor packages rsemflow needs.
#
#   Rscript scripts/install_dependencies.R          # human + mouse annotation
#   Rscript scripts/install_dependencies.R --rat    # also org.Rn.eg.db

args <- commandArgs(trailingOnly = TRUE)

cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(cran_repo) || identical(cran_repo, "@CRAN@")) {
  cran_repo <- "https://cloud.r-project.org"
}

installed <- function(pkgs) {
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
}

cran <- c("data.table", "msigdbr", "testthat")
missing_cran <- cran[!installed(cran)]
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
if ("--rat" %in% args) bioc <- c(bioc, "org.Rn.eg.db")
missing_bioc <- bioc[!installed(bioc)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

still_missing <- c(cran, bioc)[!installed(c(cran, bioc))]
if (length(still_missing)) {
  stop("Could not install: ", paste(still_missing, collapse = ", "), call. = FALSE)
}
