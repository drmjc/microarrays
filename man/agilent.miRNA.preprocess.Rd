\name{agilent.miRNA.preprocess}
\alias{agilent.miRNA.preprocess}
\title{Preprocess Agilent miRNA data.
From a table of imported Agilent data,
 add an offset, truncate, log-transform, normalize and plot a set of
Agilent microarray files. This combines \code{\link{transform_Agilent_log}},
\code{\link{agilent.miRNA.normalise.GX10}} or \code{\link[limma]{normalizeBetweenArrays}},
and \code{\link{agilent.miRNA.filter.GX10}}}
\usage{
  agilent.miRNA.preprocess(data, offset = 8,
    min.thresh = 1, percentile = 0.75, min.Pcount = 1,
    species = "mmu", plot = FALSE, do.par = TRUE,
    norm.method = c("percentile", "none", "scale", "quantile"))
}
\arguments{
  \item{data}{a \code{list} of imported Agilent data. see
  \code{\link{import.agilent}}}

  \item{offset}{add an offset. see
  \code{transform_Agilent_log}}

  \item{min.thresh}{if any values are < this value, then
  truncate them to this value.  see
  \code{\link{transform_Agilent_log}}}

  \item{percentile}{The percentile to perform
  percentile-shift normalization. only used if
  norm.method=="percentile", see
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

  \item{norm.method}{The normalisation method to use. one
  of \dQuote{percentile}, \dQuote{none}, \dQuote{scale},
  \dQuote{quantile}. If \dQuote{percentile}, then make sure
  you set the \code{\link{percentile}} parameter
  appropriately. If \dQuote{scale} or \dQuote{quantile},
  then the \code{\link[limma]{normalizeBetweenArrays}} from
  \code{limma} is used. If \dQuote{none} then no
  normalisation is performed.}
}
\value{
  a list of log-transformed, filtered, normalised Agilent
  microRNA data.
}
\description{
  Preprocess Agilent miRNA data. From a table of imported
  Agilent data, add an offset, truncate, log-transform,
  normalize and plot a set of Agilent microarray files.
  This combines \code{\link{transform_Agilent_log}},
  \code{\link{agilent.miRNA.normalise.GX10}} or
  \code{\link[limma]{normalizeBetweenArrays}}, and
  \code{\link{agilent.miRNA.filter.GX10}}
}
\author{
  Mark Cowley, 2011-08-02
}

