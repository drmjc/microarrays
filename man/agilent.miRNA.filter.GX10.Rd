\name{agilent.miRNA.filter.GX10}
\alias{agilent.miRNA.filter.GX10}
\title{Filter Agilent miRNA data
Keep only the miRNA's that start with a certain pattern, such as \dQuote{hsa}, \dQuote{mmu},
\dQuote{rno}, etc..., and only those miRNA's that are detected >= min times.}
\usage{
  agilent.miRNA.filter.GX10(data, min = 1, species = "mmu")
}
\arguments{
  \item{data}{a set of normalised Agilent miRNA data. see
  \code{\link{agilent.miRNA.normalise.GX10}}}

  \item{min}{Retain probes detected in at least \code{min}
  samples}

  \item{species}{the three letter code for which miRNA's to
  keep. Set to \code{NULL} to include all miRNA's}
}
\description{
  Filter Agilent miRNA data Keep only the miRNA's that
  start with a certain pattern, such as \dQuote{hsa},
  \dQuote{mmu}, \dQuote{rno}, etc..., and only those
  miRNA's that are detected >= min times.
}
\author{
  Mark Cowley, 2009-07-22
}

