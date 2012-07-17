\name{transform_Agilent_log}
\alias{transform_Agilent_log}
\title{Log transform and offset of Agilent data.}
\usage{
  transform_Agilent_log(data, offset = 8, minVal = 1,
    base = 2)
}
\arguments{
  \item{data}{a \code{matrix}-like object of agilent data,
  or a \code{list} of agilent data with an element called
  \sQuote{value}}

  \item{offset}{the first offset that is added to all
  values.}

  \item{minVal}{if any values are small than this minVal,
  then hard-threshold them to this minVal}

  \item{base}{the logarithm base. default=2 for log-base-2}
}
\value{
  if a matrix-like object was passed into this method, then
  a log transformed matrix is returned. If a list was
  passed into this method, then the value and error will be
  logged (base 2).
}
\description{
  Usually, >90\% of the negative values are in [-8,0), with
  the remainder being < 8. This method adds 8 to all
  values, then hard thresholds the remaining values that
  are still negative to some small positive number (eg 1)
}
\author{
  Mark Cowley, 2008-08-01
}

