\docType{methods}
\name{average.pval.replicates}
\alias{average.pval.replicates}
\alias{average.pval.replicates,ANY,missing-method}
\alias{average.pval.replicates,data.frame,ANY-method}
\alias{average.pval.replicates,matrix,character-method}
\alias{average.pval.replicates,matrix,factor-method}
\title{average.pval.replicates}
\arguments{
  \item{x}{a \code{matrix} or \code{data.frame} of p-values
  in [0,1]}

  \item{classes}{a \code{character}, \code{numeric} or
  \code{factor} vector of classes that each sample belongs
  to. in this instance, class usually represents the
  original sample from which each array has been generated
  from.}
}
\value{
  a \code{matrix} or \code{data.frame} of 'averaged'
  p.values. See Details.
}
\description{
  average.pval.replicates
}
\details{
  Given a 2D table of pvalues, such as the detection pvalue
  for each probe for each sample, compute an average
  p-value for those samples which are replicates of each
  other.\cr This uses the more correct calculation of
  \code{p*(1-ln(p)), where p = p1*p2}, which is easily
  generalised to >2 replicates.
}
\note{
  Since detection pvals=0.00000 are quite common, and
  log(0) == NaN, if \emph{any} of the p-values are 0, then
  the averaged p-value is also set to 0. I'm not 100\% sure
  this is the best approach -- it's certainly correct of
  all values are 0, but theoretically, 1 sample could have
  P=1.0 and the other P=0.
}
\author{
  Mark Cowley, Mark Pinese
}

