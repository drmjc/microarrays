\name{plot_arrayWeights}
\alias{plot_arrayWeights}
\title{Plot array weights}
\usage{
  plot_arrayWeights(aw = NULL, data = NULL, design = NULL,
    verbose = TRUE, main = "Array Weights", ...)
}
\arguments{
  \item{aw}{a named, numeric vector of array weights
  [optional]}

  \item{data}{if \code{aw=NULL}, then you can specify an
  object that must be compatible with arrayWeights. see
  \code{\link[limma]{arrayWeights}}}

  \item{design}{if \code{aw=NULL}, you can specify the
  design matrix for calculating the array weights. defaults
  to the unit vector.}

  \item{verbose}{logical}

  \item{main}{the plot title}

  \item{\dots}{arguments passed to \code{\link{barplot}}}
}
\value{
  creates a barplot of array weights, and invisibly returns
  the array weights.
}
\description{
  Function to plot (and optionally calculate) the array
  weights using \code{\link[limma]{arrayWeights}} in
  \code{limma}.
}
\author{
  Mark Cowley, 2008-07-11
}

