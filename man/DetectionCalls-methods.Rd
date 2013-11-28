\docType{methods}
\name{DetectionCalls}
\alias{DetectionCalls}
\alias{DetectionCalls,ExpressionSet,missing-method}
\alias{DetectionCalls,ExpressionSet,numeric-method}
\alias{DetectionCalls,LumiBatch,missing-method}
\alias{DetectionCalls,LumiBatch,numeric-method}
\title{DetectionCalls}
\arguments{
  \item{x}{an ExpressionSet, or LumiBatch object}

  \item{thresh}{a numeric(2) where the values represent the
  p-value thresholds for Present and Marginal,
  respectively. Some guidance re Illumina detection pvals
  here
  \url{http://www.illumina.com/Documents/products/technotes/technote_gene_expression_data_quality_control.pdf}}
}
\value{
  a character matrix, same dimension as x
}
\description{
  From microarray data, extract the detection p-values, and
  convert to P/M/A calls, ie Present/Marginal/Absent.
}
\section{Affymetrix Genome/IVT arrays}{
  Older generation Affymerix arrays use
  \code{affy::\link{mas5calls}}, where the thresholds are
  \code{c(0.04, 0.06)}.
}

\section{Affymetrix ST arrays}{
  ST arrays, that use the DABG procedure for detection
  pvalues need to have quite stringent p-value thresholds.
  We recommend \code{c(1e-05, 0.001)}.
}

\section{LumiBatch objects}{
  As far as I can tell, BeadStudio and GenomeStudio differ
  in how they report detection p-values. BeadStudio (older)
  has low values indicating expression; GenomeStudio
  reports 1-p, thus high values indicate high expression.
  This is true at least for BeadStudio 1.4.0.1, and
  GenomeStudio 1.8.0 files. In either case, the values for
  thresh should be small (default=c(0.01, 0.05)), and for
  GenomeStudio, the logic will be reversed such that probes
  with p-val > 0.95 are Marginal, and p-val > 0.99 are
  Present.
}
\examples{
require(lumi)
data(example.lumi)
head(DetectionCalls(example.lumi))
}
\author{
  Mark Cowley
}

