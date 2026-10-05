\name{compute_correlation}
\alias{compute_correlation}
\title{Compute sample correlation}
\usage{
compute_correlation(
  study,
  transform = c("vst", "log2-tpm"),
  method = c("pearson", "spearman")
)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{transform}{Expression transformation.}
\item{method}{Correlation method.}
}
\value{A sample-by-sample numeric matrix.}
\description{Computes sample correlation from transformed expression.}
