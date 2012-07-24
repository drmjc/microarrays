\name{catplot.RankProd}
\alias{catplot.RankProd}
\title{CAT plot comparing N RankProd analyses.}
\usage{
  catplot.RankProd(..., names = LETTERS)
}
\arguments{
  \item{\dots}{at least 1 rankprod object}

  \item{names}{the names of each rankprod object. length >=
  number of \dots arguments}
}
\value{
  none. a CAT plot is made
}
\description{
  All pair-wise combinations of CAT plots are made between
  each of the supplied lists produced by RP or RPAdvance.
}
\author{
  Mark Cowley, 2008-12-10
}

