\name{normalize_expression}
\alias{normalize_expression}
\title{Normalize or transform RSEMflow expression}
\usage{
normalize_expression(
  study,
  method = c("vst", "normalized-counts", "log2-tpm")
)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{method}{Transformation method.}
}
\value{A gene-by-sample numeric matrix.}
\description{Provides DESeq2 VST, normalized counts, or log2(TPM+1).}
