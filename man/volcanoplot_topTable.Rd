\name{volcanoplot_topTable}
\alias{volcanoplot_topTable}
\title{A topTable volcano plot.
Volcano plot of a toptable - similar to \code{\link[limma]{volcanoplot}} from \code{limma},
 but more flexible in terms of what is plotted on the y-axis.}
\usage{
  volcanoplot_topTable(tt, yaxis = c("p", "q", "B", "t"),
    xaxis = c("logFC", "signedFC", "FC"), yThresh = 0.05,
    lfc = log2(1.5), highlight = 0, names = tt$ID,
    main = "Volcano Plot", cex.points = 0.2, ablines = "",
    abline.col = "grey", colour.scheme = "red",
    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
    do.par = TRUE, ...)
}
\arguments{
  \item{tt}{a \code{data.frame} object from calling
  \code{\link[limma]{topTable}} with all genes.}

  \item{yaxis}{choose one of \dQuote{p}, \dQuote{q},
  \dQuote{B}, \dQuote{t}, \dQuote{absT} to plot the raw
  P-values, the q value (either from a column called
  \dQuote{q}, or \dQuote{adj.P.Val} in that order), the B
  statistics (log odds), the t-statistic, or the absolute
  t-statistic, respectively.}

  \item{xaxis}{the x-axis type. one of: \dQuote{logFC}
  (default), \dQuote{signedFC} (see
  \code{\link{logFC2signedFC}}), or \dQuote{FC} (unlogged
  Fold Change)}

  \item{yThresh}{For colouring significant points. For eg,
  choose 0.05 for P/q, or 3 for B, or 5 for t}

  \item{lfc}{the absolute log fold change threshold. eg
  0.585 for 1.5 FC, or 1.0 for a 2-fold change}

  \item{highlight}{the number of points to highlight by
  name}

  \item{names}{the names to use if highlight > 0. This
  should have as many rows as were passed into
  \code{lmFit}, and should be in the same order as the
  genes passed into lmFit. The tt's row names (which are
  numeric) will be used to index into this vector of
  names.}

  \item{main}{the plot title.}

  \item{cex.points}{see \code{\link{par}}}

  \item{ablines}{whether to add vertical dashed lines at
  \code{yThresh} and \code{xThresh}.  \dQuote{},
  \dQuote{x}, \dQuote{y}, or \dQuote{xy} for none, x (ie
  vertical lines at +/- xThresh only), y (ie horizontal
  line at yThresh only), x and y}

  \item{abline.col}{the colour of the ablines, if
  \code{ablines != ""}}

  \item{colour.scheme}{Control the higlighting of DE genes,
  via either black/red dots (\dQuote{red}) or closed/open
  circles(\dQuote{black})}

  \item{xlim}{see \code{\link{par}} if \code{NULL} it will
  be made symmetrical so that all data points fit}

  \item{ylim}{see \code{\link{par}}. if \code{NULL} it will
  fit all data points}

  \item{xlab}{The x-axis label. if \code{NULL} it's
  auto-determined from \code{xaxis}}

  \item{ylab}{The y-axis label. if \code{NULL} it's
  auto-determined from \code{yaxis}}

  \item{do.par}{logical: setup the plot parameters}

  \item{\dots}{further arguments passed to
  \code{\link{plot}}}
}
\value{
  nothing
}
\description{
  A topTable volcano plot. Volcano plot of a toptable -
  similar to \code{\link[limma]{volcanoplot}} from
  \code{limma}, but more flexible in terms of what is
  plotted on the y-axis.
}
\author{
  Mark Cowley, 2008-07-25
}

