\name{import.illumina.bgx}
\alias{import.illumina.bgx}
\title{Import an Illumina microarray manifest file.}
\usage{
  import.illumina.bgx(file, what = c("probes", "controls"))
}
\arguments{
  \item{file}{the path to the Illimina BGX file}

  \item{what}{what probes to import? one of \dQuote{probes}
  or \dQuote{controls}}
}
\value{
  a \code{data.frame} representation of the probes from a
  BGX file.
}
\description{
  Illumina microarray manifest files describe the contents
  of each microarray, including the probe names and
  sequences among many other things. You can download them
  from [1] or [2]. You should download the text versions of
  these files, unzip them and import using this function.
  You can either import the probe information,or the
  control information.
}
\note{
  This seems like a way to unzip BGX manifest files which i
  don't remember being possible:\cr if the file ends in
  zip, then unzip it. then unzip the bgx file using
  \code{gunzip -S .bgx xxxxx}, then rename it to end with
  \code{.bgx} again & import it using this function.
}
\author{
  Mark Cowley, 2008-10-23
}
\references{
  [1] \url{http://www.switchtoi.com/annotationfiles.ilmn}
  [2]
  \url{http://www.switchtoi.com/annotationprevfiles.ilmn}
}

