\docType{methods}
\name{plot_exprs_vs_detection}
\alias{plot_exprs_vs_detection}
\alias{plot_exprs_vs_detection,ExpressionSet,character,character,missing,missing,missing-method}
\alias{plot_exprs_vs_detection,ExpressionSet,character,character,numeric,character,character-method}
\alias{plot_exprs_vs_detection,ExpressionSet,character,character,numeric,character,NULL-method}
\alias{plot_exprs_vs_detection,ExpressionSet,character,character,numeric,missing,missing-method}
\alias{plot_exprs_vs_detection,ExpressionSet,character,missing,missing,missing,missing-method}
\alias{plot_exprs_vs_detection,ExpressionSet,missing,missing,missing,missing,missing-method}
\title{plot_exprs_vs_detection}
\arguments{
  \item{x}{an ExpressionSet object}

  \item{probes}{a character vector of at least 1 probe ID.
  must be in \code{featureNames(x)}}

  \item{samples}{an optional vector of sample names to
  highlight in the plot}

  \item{detection.thresh}{an optional detection threshold
  for calling 'detected' genes}

  \item{mains}{an optional plot title(s). If missing, then
  the plot title is set to gene - probe - n=<n>, where gene
  is from the \code{fData(x)[,gene.column]}, probe is from
  \code{probes}, and <n> is the number of samples in
  \code{x}.}

  \item{gene.column}{the column name within \code{fData(x)}
  which contains the gene symbols.}
}
\value{
  nothing; a 2-panel plot is produced
}
\description{
  From a set of expression data, make one 2-panel-plot for
  each probe in \code{probes}, where the top panel is the
  sorted expression levels, from low to high, and the
  bottom panel is the detection p-values. Datapoints are
  \dQuote{x} or \dQuote{o} depdending on whether they are
  above or below the \code{detection.thresh}. Optionally,
  certain named \code{samples} can the highlighted in red.
  A percentile axis is added to the top of the top panel.
  If you supply more than 1 probe, then multiple plots will
  be displayed. If you leave \code{mains=missing}, then you
  must specify \code{gene.column}, as the column name
  within \code{fData(x)[,gene.column]} which contains gene
  symbols. It's been designed to display large numbers of
  samples.
}
\examples{
if( require(lumi) ) {
  data(example.lumi)
  example.lumi
  plot_exprs_vs_detection(example.lumi, probes="oZsQEQXp9ccVIlwoQo", samples=c("A01", "A02"), 0.01, gene.column="TargetID")
}
}
\author{
  Mark Cowley
}

