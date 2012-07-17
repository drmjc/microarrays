\name{summarise.topTable}
\alias{summarise.topTable}
\title{Summarise a topTable object}
\usage{
  summarise.topTable(tt, p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    logFC.thresh = c(0.585, 1, 2))
}
\arguments{
  \item{tt}{a \code{\link[limma]{topTable}} or
  \code{topTableQ} object}

  \item{p.thresh}{a numberic vector of p-value threhshlds
  to use. set to NULL if not required.}

  \item{q.thresh}{a numeric vector of q-value threhsolds to
  use. set to NULL if not required}

  \item{logFC.thresh}{a numeric vector of absolute logFC
  threhsolds to use.  set to NULL if not required}
}
\value{
  a 1 x n \code{data.frame}, where n = np + nq where np is
  the number of p-value threhsolds, and nq is the number of
  q-value thresholds.  See also:
  \code{\link{export.DEgenes.topTable}}
}
\description{
  Get the DE genes from a \code{topTable} or
  \code{topTableQ} object, looking in the \dQuote{P.Value},
  and \dQuote{q} or \dQuote{logFC} column, returning a
  \code{data.frame} with 1 row, with numbers of DE genes in
  each column.
}
\author{
  Mark Cowley, 25/2/2008
}

