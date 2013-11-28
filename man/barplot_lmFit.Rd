\name{barplot_lmFit}
\alias{barplot_lmFit}
\title{Multipanel barplot of results from an lmFit}
\usage{
  barplot_lmFit(fit1, fit2 = NULL, data, calls = NULL,
    data.type = c("1colour", "2colour")[1], probes = NULL,
    tt = NULL, number = 10, probe2genesymbol = NULL,
    fit1.colour = "#53406A", fit2.colour = "#4F81BD",
    data.colour = NULL, hgrid.col = "black", do.par = TRUE,
    drop.fit1.intercept = FALSE,
    legend.pos = "bottomright")
}
\arguments{
  \item{fit1}{objects from lmFit. fit2 is optional.}

  \item{fit2}{objects from lmFit. fit2 is optional. fit1
  can also be just a data.frame which is useful for paired
  analyses, where you often convert the expression data
  (2n) into expression ratios (1n), prior to then fitting a
  linear model. In this instance, ags should be
  fit1=ratios, fit2=lmFitXYZ, data=rma.}

  \item{data}{data.frame of expression level data}

  \item{calls}{optional data.frame of calls, same dim as
  rma}

  \item{data.type}{\dQuote{1colour} or \dQuote{2colour}}

  \item{probes}{optional vector of probe indices, or
  probeset ids}

  \item{tt}{optional toptable of results. if supplied, you
  should set number to some}

  \item{number}{optional toptable of results. if supplied,
  you should set number to some positive integer
  corresponding to number of genes to plot.}

  \item{probe2genesymbol}{2 column table with probe ID's
  and gene symbols, respectively}

  \item{fit1.colour}{optional vector of colours for the N
  columns in fit1.  defaults to grey}

  \item{fit2.colour}{optional vector of colours for the N
  columns in fit1.  defaults to grey}

  \item{data.colour}{optional vector of colours for the N
  columns in fit1.  defaults to grey. this is ignored if
  calls != NULL}

  \item{hgrid.col}{do you want horizontal grid lines? NULL
  means no, otherwise choose a single colour.}

  \item{do.par}{logical: set the layout and the par
  settings?}

  \item{drop.fit1.intercept}{logical: drop the intercept
  term in the first fit object?}

  \item{legend.pos}{Position of the legend. See
  \code{\link{legend}}, Default = \dQuote{bottomright}}
}
\value{
  none
}
\description{
  Very flexible function to barplot results from an lmFit.
  It handles results from 2 styles of limma analysis:\cr 1:
  \dQuote{Standard analysis}: model.matrix > lmFit > eBayes
  > topTable\cr 2: \dQuote{Constrast analysis}:
  model.matrix > lmFit > fit.constrasts > eBayes >
  topTable\cr
}
\details{
  1. \dQuote{Standard analysis}\cr for each probe, do a
  barplot of the normalised data, then an errorbar plot
  utilising the coefficients and the standard errors
  (\code{stdev.unscaled * sigma}) from the lmFit1 object.
  2. \dQuote{Contrast analysis}\cr for each probe, make 3
  barplots. The first 2 are same as standard analysis, the
  3rd is an errorbar plot based on fit2 object which you
  get after doing a \code{\link[limma]{contrasts.fit}}
}
\section{Improvements}{
  Probe selection:\cr 1. probe=a numeric vector of row
  indices into the lmFit (ie same row order as data)\cr 2.
  probe=vector of probesetID's which are in the rownames of
  data and fit1 [and fit2] [and calls]\cr 3. supply a
  topTable object, and set the number of rows from top to
  bottom to plot.\cr this can be from an F-test or t-test
}

\section{Colouring array data}{
  if you supply a 'calls' object which is same dim as data,
  and contains \dQuote{P}, \dQuote{M} or \dQuote{A}, then
  the bars for the expression data will be coloured green,
  orange or red, respectively.
}
\examples{
\dontrun{
barplot_lmFit(fit1, data=rma, probes=c("10543233", "10411107"))
barplot_lmFit(fit1, fit2, data=rma, probes=c("10543233", "10411107"))
barplot_lmFit(fit1, data=rma, tt=topTable, number=2)
barplot_lmFit(fit1, fit2, data=rma, tt=topTable, number=2)
barplot_lmFit(fit1, fit2, data=rma, calls=calls, tt=topTable, number=2)
}
}
\author{
  Mark Cowley, 2009-07-16
}

