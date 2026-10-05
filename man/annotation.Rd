\name{annotate_ensembl}
\alias{annotate_ensembl}
\alias{annotate_gtf}
\title{Annotate RSEMflow genes}
\usage{
annotate_ensembl(study, species = "auto")

annotate_gtf(study, gtf)
}
\arguments{
\item{study}{An \code{rsemflow_study} or study path.}
\item{species}{\code{auto}, \code{human}, \code{mouse}, or \code{rat}.}
\item{gtf}{Path to a GTF or GTF.GZ file.}
}
\value{A gene annotation data frame retaining the original gene ID.}
\description{Map versioned Ensembl IDs to stable IDs and gene annotation.}
