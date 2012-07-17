\name{summarise.topTable.updown}
\alias{summarise.topTable.updown}
\title{Summarise a toptable into the up/down regulated genes}
\usage{
  summarise.topTable.updown(tt,
    p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    logFC.thresh = c(0.585, 1, 2))
}
\arguments{
  \item{tt}{a limma \code{topTable}, or
  \code{list(topTable)}}

  \item{p.thresh}{a vector of p value thresholds, set to
  \code{NULL} to ignore}

  \item{q.thresh}{a vector of adj.P.Val thresholds, set to
  \code{NULL} to ignore.}

  \item{logFC.thresh}{a vector of absolute logFC value
  thresholds}
}
\value{
  a data.frame of counts of differentially expressed genes
}
\description{
  Summarise a toptable into the up/down regulated genes
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{summarise.topTable}},
  \code{\link{summarise.topTable.list}}
}

