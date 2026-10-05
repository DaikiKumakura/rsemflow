\name{compute_pca}
\alias{compute_pca}
\title{Compute PCA from RSEMflow expression}
\usage{
compute_pca(
  study,
  transform = c("vst", "log2-tpm"),
  ntop = 500L,
  components = 10L
)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{transform}{Expression transformation.}
\item{ntop}{Number of most variable genes.}
\item{components}{Maximum number of PCs returned.}
}
\value{A list with \code{scores}, \code{loadings}, and \code{variance}.}
\description{Computes tabular PCA outputs without visualization.}
