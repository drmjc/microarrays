\name{summarise.topTable.list}
\alias{summarise.topTable.list}
\title{Summarise a list of topTables.}
\usage{
  summarise.topTable.list(tt.list,
    p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    logFC.thresh = c(0.585, 1, 2))
}
\arguments{
  \item{tt.list}{must be a list of toptable objects, and no
  F stat tables. all need to have logFC columns.}

  \item{p.thresh}{a numberic vector of p-value threhshlds
  to use. set to NULL if not required.}

  \item{q.thresh}{a numeric vector of q-value threhsolds to
  use. set to NULL if not required}

  \item{logFC.thresh}{a numeric vector of absolute logFC
  threhsolds to use.  set to NULL if not required}
}
\value{
  a \code{data.frame} See also:
  \code{\link{summarise.topTable}},
  \code{\link{summarise.topTable.updown}}
}
\description{
  This represnts a new style of summary, where N top tables
  makes N columns (plus some descriptive ones...)
}
\author{
  Mark Cowley, 2009-12-21
}

