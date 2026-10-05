\name{read_study}
\alias{read_study}
\alias{write_study}
\title{Read and write RSEMflow studies}
\usage{
read_study(path)

write_study(study, path)
}
\arguments{
\item{path}{Study directory or \code{study.rds} for reading; output directory for writing.}
\item{study}{An \code{rsemflow_study}.}
}
\value{\code{read_study()} returns an \code{rsemflow_study}; \code{write_study()} returns the output path invisibly.}
\description{Serialize or restore a local RSEMflow study.}
