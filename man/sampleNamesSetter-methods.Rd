\name{sampleNames<-}
\alias{sampleNames<-}
\alias{sampleNames<-,LumiBatch-method}
\title{set sampleNames in LumiBatch objects}
\arguments{
  \item{object}{an \code{LumiBatch} object}

  \item{value}{a character vector of column names}
}
\value{
  a \code{LumiBatch} object
}
\description{
  set sampleNames in LumiBatch objects, including the QC,
  vstParameter, transformFun and controlData slots
}
\note{
  code came from lumi_2.8.0 & modified by MJC to set the
  controlData slot properly
}

