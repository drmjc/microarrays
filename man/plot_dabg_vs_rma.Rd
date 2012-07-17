\name{plot_dabg_vs_rma}
\alias{plot_dabg_vs_rma}
\title{Density plot of DABG vs expression levels
Plot the RMA expression levels as density plots, split by various DABG
P-value thresholds.}
\usage{
  plot_dabg_vs_rma(dabg, rma,
    thresholds = c(0.05, 0.01, 0.001, 1e-04, 1e-05))
}
\arguments{
  \item{dabg}{a \code{data.frame} of DABG P-Values}

  \item{rma}{a \code{data.frame} of expression levels
  (usually RMA normalised levels, thus the name)}

  \item{thresholds}{a vector of DABG p-value thresholds}
}
\description{
  Density plot of DABG vs expression levels Plot the RMA
  expression levels as density plots, split by various DABG
  P-value thresholds.
}
\examples{
\dontrun{
plot_dabg_vs_rma(dabg, rma)
}
}
\author{
  Mark Cowley, 17/1/08
}

