\name{subset.eSet}
\alias{subset.eSet}
\title{subset an eSet object
subset the rows (features) and columns(samples) within an eSet object.}
\usage{
  \method{subset}{eSet} (x, subset, select, drop = FALSE,
    ...)
}
\arguments{
  \item{x}{an eSet}

  \item{subset}{logical expression indicating elements or
  rows to keep: missing values are taken as \code{FALSE}.}

  \item{select}{logical expression, indicating columns to
  select from a \code{data.frame}.}

  \item{drop}{logical: ignored}

  \item{\dots}{ignored}
}
\value{
  an eSet with fewer rows and or columns
}
\description{
  subset an eSet object subset the rows (features) and
  columns(samples) within an eSet object.
}
\note{
  if you have a subclass of eSet, there may be slots that
  are not subsetted. For example LumiBatch objets can have
  an extra controlData slot.
}
\author{
  Mark Cowley, 2011-09-01
}
\seealso{
  \code{\link{subset.LumiBatch}}
}

