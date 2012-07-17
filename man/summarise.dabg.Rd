\name{summarise.dabg}
\alias{summarise.dabg}
\title{Sumarise an Affymetrix DABG file
The Detected Above Background method by Affymetrix generates a P-Value for the probability
of each probeset being expressed above the matched control probes.
Import an Affymetrix DABG results table, and: boxplot of these p-values for each array
summary table of how many probesets pass certain P-value thresholds. see
also \code{\link{plot_dabg_vs_rma}}}
\usage{
  summarise.dabg(file)
}
\arguments{
  \item{file}{the path to an Affymetrix DABG result file.}
}
\description{
  Sumarise an Affymetrix DABG file The Detected Above
  Background method by Affymetrix generates a P-Value for
  the probability of each probeset being expressed above
  the matched control probes. Import an Affymetrix DABG
  results table, and: boxplot of these p-values for each
  array summary table of how many probesets pass certain
  P-value thresholds. see also
  \code{\link{plot_dabg_vs_rma}}
}
\examples{
\dontrun{
dabg <- summarise.dabg("./dabg.summary.txt")
}
}
\author{
  Mark Cowley, 17/1/08
}

