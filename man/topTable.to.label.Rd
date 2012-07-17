\name{topTable.to.label}
\alias{topTable.to.label}
\title{Take a toptable (Fstat or t-stat) and make a label for each row.}
\usage{
  topTable.to.label(tt, probes = NULL,
    probe2genesymbol = NULL)
}
\arguments{
  \item{tt}{a toptable}

  \item{probes}{if \code{NULL}, all probes from
  \code{tt$ID} are used. otherwise, probes can be a numeric
  vector (eg 1 10) or character vector corresponding to
  \code{tt$ID}}

  \item{probe2genesymbol}{a 2 column \code{data.frame} with
  probe ID's in column 1, and gene symbols in column 2. If
  there's no gene symbol for a probe, then there should be
  an \code{NA}.}
}
\value{
  a \code{character vector} with same length as probes, eg:
  [1] "probeset:10543233 F=1964 P=7.17e-13 FDR=1.16e-10"\cr
  [2] "probeset: 10410211 t=9.26 P=8.49e-06 FDR=0.0193"\cr
}
\description{
  Take a toptable (Fstat or t-stat) and make a label for
  each row.
}
\author{
  Mark Cowley, 2009-07-16
}

