\name{LumiBatch2gct}
\alias{LumiBatch2gct}
\title{convert a LumiBatch to a GCT object}
\usage{
  LumiBatch2gct(x, description.column = "Description")
}
\arguments{
  \item{x}{a LumiBatch object}

  \item{description.column}{the description column within
  \code{fData(x)} to grab the data for the resultant
  Description column in the GCT.}
}
\value{
  a GCT object. this is currently just a \code{data.frame}
  representation.
}
\description{
  convert a LumiBatch to a GenePattern GCT object. There
  must be an \code{fData(x)$Description} column.
}
\section{TODO}{
  - have a \code{GCT} class, and an \code{as.gct} method to
  convert the various types to \code{GCT-class} -
  LumiBatch2res, which would need a LumiBatch2calls method
}
\author{
  Mark Cowley, 2012-05-02
}

