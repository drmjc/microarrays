\name{normalizeQuantiles.list}
\alias{normalizeQuantiles.list}
\title{Function to qnorm a list of matrices It joins the list elements into one
large data.frame, then quantile normalizes this large data.frame. This then
gets deconstructed back into a list with the same properties as the supplied
list}
\usage{
  normalizeQuantiles.list(x)
}
\arguments{
  \item{x}{a \code{list} which can be converted to a
  \code{data.frame}}
}
\value{
  a \code{data.frame} of quantile normalised data
}
\description{
  Function to qnorm a list of matrices It joins the list
  elements into one large data.frame, then quantile
  normalizes this large data.frame. This then gets
  deconstructed back into a list with the same properties
  as the supplied list
}
\author{
  Mark Cowley, 23 Nov 2005
}

