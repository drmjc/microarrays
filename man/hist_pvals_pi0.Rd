\name{hist_pvals_pi0}
\alias{hist_pvals_pi0}
\title{histrogram of p-values}
\usage{
  hist_pvals_pi0(pvals, bins = 20, main = "",
    xlab = "unadjusted P-values", annotate = TRUE,
    theme = c("red", "black"))
}
\arguments{
  \item{pvals}{a numeric vector of unadjusted pvalues, in
  [0,1]}

  \item{bins}{the number of histogram bins}

  \item{main}{the plot title}

  \item{xlab}{the plot x-axis label}

  \item{annotate}{logical: if \code{TRUE}, annotate the
  plots with 2 horizontal lines, at pi0 under H0, and
  estimated pi0}

  \item{theme}{one of \dQuote{red} or \dQuote{black}, where
  \dQuote{red} colours the bars from deep red to white, or
  \dQuote{black} leaves the bars black & white.}
}
\value{
  nothing.
}
\description{
  Plot a histogram of unadjusted P-values, and the expected
  proportion of truly nulls, and the estimated pi0 from
  fitting the non-\code{NA} pvalues through qvalue with
  default options.
}
\author{
  Mark Cowley, 2008-10-02
}

