\name{read_rsem_study}
\alias{read_rsem_study}
\title{Import RSEM gene-level results}
\usage{
read_rsem_study(rsem_dir, metadata, categorical = character(), numeric = character())
}
\arguments{
\item{rsem_dir}{Directory containing \code{*.genes.results}.}
\item{metadata}{TSV or CSV metadata containing \code{sample_id}.}
\item{categorical}{Metadata columns forced to categorical.}
\item{numeric}{Metadata columns forced to numeric.}
}
\value{An object of class \code{rsemflow_study}.}
\description{Validates and imports RSEM gene-level expression with tximport.}
