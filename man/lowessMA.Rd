\name{lowessMA}
\alias{lowessMA}
\title{add a lowess curve to a plot created by plotMA}
\usage{
  lowessMA(MA, array = 1, col = "purple", ...)
}
\arguments{
  \item{MA}{a \code{list}, \code{MAList}, \code{RGList},
  \code{MArrayLM}}

  \item{array}{an index into the MA object}

  \item{col}{colur}

  \item{\dots}{arguments passed to
  \code{\link[stats]{lowess}}}
}
\description{
  add a lowess curve to a plot created by plotMA
}
\examples{
\dontrun{
plotMA(RG, 1)
lowessMA(RG, 1, col="red")
}
}
\author{
  Mark Cowley, 2008-08-19
}

