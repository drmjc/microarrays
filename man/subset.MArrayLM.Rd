\name{subset.MArrayLM}
\alias{subset.MArrayLM}
\title{Subset an MArrayLM object.
S3 method to subset the rows of an MArrayLM object, a la \code{\link[base]{subset}}}
\usage{
  subset.MArrayLM(x, subset, select, drop = FALSE, ...)
}
\arguments{
  \item{x}{an MArrayLM object to be subsetted.}

  \item{subset}{a logical vector of rows/observations to
  include. Different to the default subset, you can also
  specify a character vector, or a numeric vector of
  indices. Leave missing if you want to ignore this
  argument.}

  \item{select}{a logical vector of columns/conditions to
  include. Like subset, you can specify a character or
  numeric argument. Don't try to set this if ncol(x) == 1l
  you'll get a warning and this arg will be ignored.}

  \item{drop}{ignored.}

  \item{\dots}{ignored.}
}
\value{
  an MArrayLM with less rows, and or columns.
}
\description{
  Subset an MArrayLM object. S3 method to subset the rows
  of an MArrayLM object, a la \code{\link[base]{subset}}
}
\author{
  Mark Cowley, 2009-11-10
}
\seealso{
  \code{\link{subset}} \code{\link{subset.data.frame}}
}

