\docType{methods}
\name{ExtractionSoftware}
\alias{ExtractionSoftware}
\alias{ExtractionSoftware,LumiBatch-method}
\title{ExtractionSoftware}
\arguments{
  \item{x}{a LumiBatch object}
}
\value{
  \sQuote{GenomeStudio}, or \sQuote{BeadStudio}
}
\description{
  Determine which extraction software was used to create
  the Illumina TXT file. The first header of the TXT file
  looks something like this: \sQuote{Illumina Inc.
  GenomeStudio version 1.8.0}, or \sQuote{Illumina Inc.
  BeadStudio version 1.4.0.1}. This information is embedded
  within the LumiBatch object.
}
\examples{
require(lumi)
data(example.lumi)
ExtractionSoftware(example.lumi)
}
\author{
  Mark Cowley
}

