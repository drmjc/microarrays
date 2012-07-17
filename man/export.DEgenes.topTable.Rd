\name{export.DEgenes.topTable}
\alias{export.DEgenes.topTable}
\title{Export the DE genes from a topTable.
Get the DE genes passing various P or Q value threhsolds, and export them to
a xls file.}
\usage{
  export.DEgenes.topTable(tt, file = NULL,
    p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    lfc.thresh = c(0.585, 1, 2), values = FALSE)
}
\arguments{
  \item{tt}{a top.table, or a list of top tables}

  \item{file}{the path to the output tsv file}

  \item{p.thresh}{vectors of thresholds for p.values}

  \item{q.thresh}{vectors of thresholds for FDRs}

  \item{lfc.thresh}{vector of POSITIVE log base 2 FC
  thresholds. Note that absolute logFC's are used}

  \item{values}{also return the
  \sQuote{p}/\sQuote{q}/\sQuote{logFC} value?}
}
\description{
  Export the DE genes from a topTable. Get the DE genes
  passing various P or Q value threhsolds, and export them
  to a xls file.
}
\author{
  Mark Cowley, 25/2/2008
}
\seealso{
  \code{\link{DEgenes.topTable}}
}

