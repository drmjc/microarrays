\name{topTable.all}
\alias{topTable.all}
\title{Run topTable on all contrasts from a linear model fit
Function to create a topTable for every contrast from an lmFit object,
including the F-stat, if >1 contrasts}
\usage{
  topTable.all(fit, number = Inf, genelist = fit$genes,
    adjust.method = "BH", sort.by = "p", resort.by = NULL,
    p.value = 1, lfc = 0)
}
\arguments{
  \item{fit}{an \code{lmFit} object}

  \item{number}{the number of DE genes. default=Inf which
  is all genes tested}

  \item{genelist}{vector of genes to do the topTable on}

  \item{adjust.method}{the multiple test correction method.
  default=\dQuote{BH}.  see \code{\link[limma]{topTable}}}

  \item{sort.by}{which column to sort by? see
  \code{\link[limma]{topTable}}}

  \item{resort.by}{character string specifying statistic to
  sort the selected genes by in the output
  \code{data.frame}.  Possibilities are the same as for
  \code{sort.by}.}

  \item{p.value}{a p.value threshold. default=1.0, ie all
  genes}

  \item{lfc}{a logFC threhsold to use. default=0.0, ie all
  genes}
}
\value{
  a list of topTable objects, named by the coefficient
  names
}
\description{
  Run topTable on all contrasts from a linear model fit
  Function to create a topTable for every contrast from an
  lmFit object, including the F-stat, if >1 contrasts
}
\author{
  Mark Cowley, 2010-10-07
}

