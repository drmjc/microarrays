\name{make.sample.info.Agilent}
\alias{make.sample.info.Agilent}
\title{Autocreate sample.info from Agilent filenames
From a vector of Agilent GeneView file names, create a generic \code{sample.info}
 object.}
\usage{
  make.sample.info.Agilent(filenames, outfile = NULL)
}
\arguments{
  \item{filenames}{a character vector of filenames}

  \item{outfile}{if \code{NULL} or \code{FALSE}, then do
  not write this to a file. if \code{TRUE}, then write to
  \dQuote{./sample.info.xls}, otherwise write to the
  specified file.}
}
\value{
  a \code{data.frame} with these columns:
  \dQuote{FileName}, \dQuote{Barcode}, \dQuote{SlideID},
  \dQuote{ArrayID}, \dQuote{SampleID}
}
\description{
  Autocreate sample.info from Agilent filenames From a
  vector of Agilent GeneView file names, create a generic
  \code{sample.info} object.
}
\author{
  Mark Cowley, 2008-08-01
}

