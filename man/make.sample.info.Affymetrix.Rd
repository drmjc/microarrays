\name{make.sample.info.Affymetrix}
\alias{make.sample.info.Affymetrix}
\title{From a vector of CEL file names, create a generic sample.info object.}
\usage{
  make.sample.info.Affymetrix(files = NULL, outfile = NULL,
    auto.batch = TRUE, batch.thresh = 1440)
}
\arguments{
  \item{files}{a character vector of filenames}

  \item{outfile}{if NULL or FALSE, then do write this to a
  file if TRUE, then write to ./sample.info.xls, otherwise
  write to the specified file.}

  \item{auto.batch}{try to determine any processing
  batches. see batch.thresh}

  \item{batch.thresh}{if auto.batch, then increment the
  numeric batch ID if the time between two consecutive
  files is greater than this threshold in MINTUES. Default
  is 1 day/}
}
\value{
  a data.frame
}
\description{
  From a vector of CEL file names, create a generic
  sample.info object.
}
\author{
  Mark Cowley, 2008-07-29
}

