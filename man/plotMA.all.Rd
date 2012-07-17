\name{plotMA.all}
\alias{plotMA.all}
\title{Make N MAplots, writing out files into a directory.}
\usage{
  plotMA.all(MA, dir, width = 1024, height = 768, ...)
}
\arguments{
  \item{MA}{an MAList object}

  \item{dir}{the output directory}

  \item{width}{arguments to png}

  \item{height}{arguments to png}

  \item{\dots}{further args passed to plotMA.}
}
\value{
  none. creates N files, one per array.  See also:
  limma:plotMA
}
\description{
  Make N MAplots, writing out files into a directory.
}
\author{
  Mark Cowley, 2009-11-02
}

