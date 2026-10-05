\name{run_gsea}
\alias{run_gsea}
\title{Run preranked GSEA}
\usage{
run_gsea(
  differential,
  species = "auto",
  collection = "hallmark",
  subcollection = NULL,
  genesets = NULL,
  rank = "stat",
  min_size = 15L,
  max_size = 500L,
  eps = 0
)
}
\arguments{
\item{differential}{Differential-expression data frame or TSV path.}
\item{species}{Species name/alias.}
\item{collection}{MSigDB collection alias or collection code.}
\item{subcollection}{Optional MSigDB subcollection.}
\item{genesets}{Optional named list or GMT file.}
\item{rank}{Ranking-statistic column.}
\item{min_size}{Minimum gene-set size.}
\item{max_size}{Maximum gene-set size.}
\item{eps}{fgsea multilevel epsilon.}
}
\value{A GSEA result data frame.}
\description{Runs preranked GSEA using Ensembl identifiers and fgsea.}
