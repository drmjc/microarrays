\name{convert.partek2topTable}
\alias{convert.partek2topTable}
\title{Convert Partek to topTable
Convert a toptable made by Partek's ANOVA into a limma-like-topTable.
 I expect there
to be a few columns of annotation, including the Probeset.ID; then there
should be at least a p value, t stat, fold change column. If there's a q
value column then great it will be used, otherwise q-values will be
generated for you. Do NOT let there be > 1 pval/qval/t-stat column...}
\usage{
  convert.partek2topTable(partek)
}
\arguments{
  \item{partek}{a \code{data.frame} created by Partek, and
  imported into \R}
}
\value{
  a \code{data.frame}
}
\description{
  Convert Partek to topTable Convert a toptable made by
  Partek's ANOVA into a limma-like-topTable.  I expect
  there to be a few columns of annotation, including the
  Probeset.ID; then there should be at least a p value, t
  stat, fold change column. If there's a q value column
  then great it will be used, otherwise q-values will be
  generated for you. Do NOT let there be > 1
  pval/qval/t-stat column...
}
\author{
  Mark Cowley, 2009-12-21
}

