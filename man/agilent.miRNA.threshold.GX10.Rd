\name{agilent.miRNA.threshold.GX10}
\alias{agilent.miRNA.threshold.GX10}
\title{Threshold Agilent miRNA data
Threshold raw Agilent miRNA data such that all values < thresh become =
thresh Optionally also log-base-2 transform the results afterwards}
\usage{
  agilent.miRNA.threshold.GX10(data, thresh = 1,
    log = TRUE)
}
\arguments{
  \item{data}{a \code{list} of imported Agilent data. see
  \code{\link{import.agilent}}}

  \item{thresh}{if any values are < this value, then
  truncate them to this value.}

  \item{log}{logical: perform log-base-2 transformation on
  the truncated data?}
}
\description{
  Threshold Agilent miRNA data Threshold raw Agilent miRNA
  data such that all values < thresh become = thresh
  Optionally also log-base-2 transform the results
  afterwards
}
\author{
  Mark Cowley, 2009-07-22
}

