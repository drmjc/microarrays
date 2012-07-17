\name{RankProd.signed.score}
\alias{RankProd.signed.score}
\title{From the results of a RankProd analysis, create a signed score suitable for
something like GSEA analysis.}
\usage{
  RankProd.signed.score(rp,
    method = c("inv.log", "inverse"))
}
\arguments{
  \item{rp}{a rankprod result, or a list containing a table
  called \dQuote{RPs}}

  \item{method}{one of \dQuote{inv.log}, or
  \dQuote{inverse}}
}
\value{
  a named vector of rank product signed scores
}
\description{
  From the results of a RankProd analysis, create a signed
  score suitable for something like GSEA analysis.
}
\author{
  Mark Cowley, 2009-01-09
}

