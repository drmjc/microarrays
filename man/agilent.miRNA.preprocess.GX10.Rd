\name{agilent.miRNA.preprocess.GX10}
\alias{agilent.miRNA.preprocess.GX10}
\title{Preprocess Agilent miRNA data
Preprocess raw Agilent data, then threshold it to a minimum value, log
transform, normalise each array such that the 75th percentile is 0.0, throw
away miRs from other species, and keep those that are detected at least
min.Pcount times. Optionally makes a 2x2 portrait layout plot @ each step.
This is the default procedure by GeneSpring GX 10 (not sure about newer versions)}
\usage{
  agilent.miRNA.preprocess.GX10(data, min.thresh = 1,
    percentile = 0.75, min.Pcount = 1, species = "mmu",
    plot = FALSE, do.par = TRUE)
}
\arguments{
  \item{data}{a \code{list} of imported Agilent data. see
  \code{\link{import.agilent}}}

  \item{min.thresh}{if any values are < this value, then
  truncate them to this value.  see
  \code{\link{agilent.miRNA.threshold.GX10}}}

  \item{percentile}{The percentile to perform
  percentile-shift normalization. See
  \code{\link{agilent.miRNA.normalise.GX10}}}

  \item{min.Pcount}{integer exclude genes detected in fewer
  than \code{min.Pcount} samples. see
  \code{\link{agilent.miRNA.filter.GX10}}}

  \item{species}{the 3 letter code for which species of
  miRNA's to include.  see
  \code{\link{agilent.miRNA.filter.GX10}}}

  \item{plot}{logical: add boxplots of data during the
  transformation process}

  \item{do.par}{logical: configure the plotting device
  settings?}
}
\value{
  a set of thresholded, log-transformed, filtered,
  normalized Agilent miRNA data
}
\description{
  Preprocess Agilent miRNA data Preprocess raw Agilent
  data, then threshold it to a minimum value, log
  transform, normalise each array such that the 75th
  percentile is 0.0, throw away miRs from other species,
  and keep those that are detected at least min.Pcount
  times. Optionally makes a 2x2 portrait layout plot @ each
  step. This is the default procedure by GeneSpring GX 10
  (not sure about newer versions)
}
\author{
  Mark Cowley, 2009-07-22
}

