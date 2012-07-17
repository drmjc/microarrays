\name{DEgenes.topTable}
\alias{DEgenes.topTable}
\title{Get the DE genes from a topTable.
Get the lists of DE genes that pass various P or Q thresholds, see also
\code{\link{summarise.topTable}}}
\usage{
  DEgenes.topTable(tt, p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    lfc.thresh = c(0.585, 1, 2), values = FALSE,
    as.df = FALSE, direction = c("either", "up", "down"))
}
\arguments{
  \item{tt}{a top.table, or a list of top tables}

  \item{p.thresh}{vectors of thresholds for p.values}

  \item{q.thresh}{vectors of thresholds for FDRs}

  \item{lfc.thresh}{vector of POSITIVE log base 2 FC
  thresholds. Note that absolute logFC's are used}

  \item{values}{also return the
  \sQuote{p}/\sQuote{q}/\sQuote{logFC} value?}

  \item{as.df}{logical if \code{FALSE} the result will be a
  list, otherwise will be a \code{data.frame}}

  \item{direction}{Which direction can the gene change be?
  one of \dQuote{either}, \dQuote{up} or \dQuote{down}}
}
\value{
  A list of DE genes from a toptable
}
\description{
  Get the DE genes from a topTable. Get the lists of DE
  genes that pass various P or Q thresholds, see also
  \code{\link{summarise.topTable}}
}
\author{
  Mark Cowley, 25/2/2008
}
\seealso{
  \code{\link{summarise.topTable}}
}

