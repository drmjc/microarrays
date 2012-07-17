\name{volcanoplot_topTable_PandQ}
\alias{volcanoplot_topTable_PandQ}
\title{A wrapper function to plot 2 volcano plots side by side, for the logFC vs P
and q values.}
\usage{
  volcanoplot_topTable_PandQ(tt, pThresh = 0.001,
    qThresh = 0.05, lfc = log2(1.5), cex.points = 0.5,
    ablines = "", ...)
}
\arguments{
  \item{tt}{a \code{data.frame} object from calling
  \code{\link[limma]{topTable}} with all genes.}

  \item{pThresh}{For colouring significant points in the
  p-value panel., default=0.001}

  \item{qThresh}{For colouring significant points in the
  FDR panel., default=0.05}

  \item{lfc}{the absolute log fold change threshold. eg
  0.585 for 1.5 FC, or 1.0 for a 2-fold change}

  \item{cex.points}{see \code{\link{par}}}

  \item{ablines}{whether to add vertical dashed lines at
  \code{yThresh} and \code{xThresh}.  \dQuote{},
  \dQuote{x}, \dQuote{y}, or \dQuote{xy} for none, x (ie
  vertical lines at +/- xThresh only), y (ie horizontal
  line at yThresh only), x and y}

  \item{\dots}{further arguments passed to
  \code{\link{volcanoplot_topTable}}}
}
\description{
  A wrapper function to plot 2 volcano plots side by side,
  for the logFC vs P and q values.
}
\author{
  Mark Cowley, 2008-08-01
}
\seealso{
  \code{\link{volcanoplot_topTable}}
}

