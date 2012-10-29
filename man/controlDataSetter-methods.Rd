\docType{methods}
\name{controlData<-}
\alias{controlData<-}
\alias{controlData<-,LumiBatch,data.frame-method}
\alias{controlData<-,LumiBatch,NULL-method}
\title{replace the controlData slot in a LumiBatch object}
\arguments{
  \item{object}{a LumiBatch object}

  \item{value}{a a \code{data.frame} with first two columns
  as \dQuote{controlType} and \dQuote{ProbeID}. The rest
  columns are the control probe expression amplitudes for
  individual samples. To delete the slot, set
  \code{value=NULL}.}
}
\value{
  something
}
\description{
  replace the controlData slot in a LumiBatch object
}
\author{
  Mark Cowley
}

