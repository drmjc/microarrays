\name{subset.LumiBatch}
\alias{subset.LumiBatch}
\title{subset a LumiBatch object}
\usage{
  \method{subset}{LumiBatch} (x, subset, select,
    drop = FALSE, ...)
}
\arguments{
  \item{x}{an LumiBatch}

  \item{subset}{logical expression indicating elements or
  rows to keep: missing values are taken as \code{FALSE}.}

  \item{select}{logical expression, indicating columns to
  select from a \code{data.frame}.}

  \item{drop}{logical: ignored}

  \item{\dots}{ignored}
}
\value{
  an LumiBatch with fewer rows and or columns
}
\description{
  This is preferable to subset.eSet, since there are extra
  slots in a LumiBatch
}
\details{
  The support for subsetting the controlData is still
  inadequately tested. For instance, do all controlData
  slots have the first 2 columns being "?" and "ProbeID"?
}
\author{
  Mark Cowley, 2011-09-01
}

