\name{differential_deseq2}
\alias{differential_deseq2}
\title{Metadata-driven DESeq2 differential expression}
\usage{
differential_deseq2(
  study,
  formula,
  term,
  metadata = NULL,
  reference = character(),
  contrast = character(),
  at = character(),
  test = c("wald", "lrt"),
  reduced = NULL,
  categorical = character(),
  numeric = character(),
  alpha = 0.05
)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{formula}{R formula or formula string.}
\item{term}{Metadata variable or interaction term to report.}
\item{metadata}{Optional replacement metadata path or data frame.}
\item{reference}{Repeatable reference specifications such as \code{condition::Control}.}
\item{contrast}{Repeatable comparisons such as \code{condition::Drug::Control}.}
\item{at}{Optional simple-effect condition such as \code{genotype::KO}.}
\item{test}{\code{wald} or \code{lrt}.}
\item{reduced}{Reduced formula for an LRT.}
\item{categorical}{Columns forced to categorical.}
\item{numeric}{Columns forced to numeric.}
\item{alpha}{DESeq2 alpha threshold.}
}
\value{A list containing the DESeq2 object, design matrix, manifest, and result tables.}
\description{Runs DESeq2 from an explicit metadata formula with rank validation and reproducible contrasts.}
