\name{run_gsva}
\alias{run_gsva}
\title{Compute GSVA or ssGSEA pathway scores}
\usage{
run_gsva(
  study,
  species = "auto",
  collection = "hallmark",
  subcollection = NULL,
  genesets = NULL,
  method = c("gsva", "ssgsea"),
  min_size = 10L,
  max_size = 500L
)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{species}{Species name/alias.}
\item{collection}{MSigDB collection alias or collection code.}
\item{subcollection}{Optional MSigDB subcollection.}
\item{genesets}{Optional named list or GMT file.}
\item{method}{\code{gsva} or \code{ssgsea}.}
\item{min_size}{Minimum mapped gene-set size.}
\item{max_size}{Maximum mapped gene-set size.}
}
\value{A pathway-by-sample data frame.}
\description{Computes sample-level pathway scores from log2(TPM+1) using the current GSVA parameter-object API.}
