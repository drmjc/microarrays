\name{agilent.miRNA.normalise.GX10}
\alias{agilent.miRNA.normalise.GX10}
\title{Normalise Agilent miRNA data
Normalise an Agilent miRNA data set using the default method from Agilent
GeneSpring GX 10 which is percentile shift normalisation. From the GX 10
manual:\cr
18.2.1 Percentile Shift Normalization\cr
Percentile shift normalization is a global normalization, where the location
of all the spot intensities in an array are adjusted. This normalization
takes each column in an experiment independently, and computes the
percentile of the expression values for this array, across all spots (where
n has a range from 0-100 and n=50 is the median). It subtracts this value
from the expression value of each entity.}
\usage{
  agilent.miRNA.normalise.GX10(data, pctile = 0.75)
}
\arguments{
  \item{data}{a \code{list} of imported, thresholded,
  log-base-2 Agilent data.  see
  \code{\link{agilent.miRNA.threshold.GX10}}}

  \item{pctile}{the percentile to normalize to.
  default=0.75, ie the 75th percentile. Must be in [0,1]}
}
\value{
  a list of normalised miRNA data
}
\description{
  Normalise Agilent miRNA data Normalise an Agilent miRNA
  data set using the default method from Agilent GeneSpring
  GX 10 which is percentile shift normalisation. From the
  GX 10 manual:\cr 18.2.1 Percentile Shift Normalization\cr
  Percentile shift normalization is a global normalization,
  where the location of all the spot intensities in an
  array are adjusted. This normalization takes each column
  in an experiment independently, and computes the
  percentile of the expression values for this array,
  across all spots (where n has a range from 0-100 and n=50
  is the median). It subtracts this value from the
  expression value of each entity.
}
\author{
  Mark Cowley, 2009-07-22
}

