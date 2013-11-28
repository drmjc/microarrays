\name{dotplot}
\alias{dotplot}
\alias{dotplot,ExpressionSet,ExpressionSet,ANY-method}
\alias{dotplot,ExpressionSet,missing,ANY-method}
\alias{dotplot,ExpressionSet,numeric,ANY-method}
\alias{dotplot,numeric,ExpressionSet,ANY-method}
\title{dotplots on ExpressionSet objects}
\arguments{
  \item{x}{an ExpressionSet}

  \item{y}{an ExpressionSet, or missing.}

  \item{sort}{logical: if \code{TRUE}, then sort from low
  to high (default=TRUE). only applies if y is missing.}

  \item{feature}{the index or name of the feature to plot.
  default=1. In the 2D case, if the feature exists in both
  ExpressionSet's then that will be seleected. If plotting
  different features from x and y, then leave feature=1,
  and select the feature of interest from each expression
  set using '['. see examples.}

  \item{xlab}{see par}

  \item{ylab}{see par}

  \item{main}{see par}

  \item{add.mean}{logical: add a horizontal line about the
  mean of that feature}

  \item{samples}{an optional character vector of sample
  names to highlight and label in red.}

  \item{lowess.col}{the colour of the 2D loess line.
  \code{NA} to not show one, otherwise choose a named
  colour.}
}
\value{
  nothing
}
\description{
  Plot either a single feature from an ExpressionSet as a
  dotplot (default is sorted low to high), or generate a 2D
  dotplot from combinations of 2 ExpressionSet's, or
  ExpressionSet and a named numeric, where those names
  overlap with the featureNames from the ExpressionSet. As
  always, the 2D objects need to have some overlap in their
  names, but otherwise, data unique to each set will be
  silently ignored.
}
\section{2D dotplots}{
  these are just xy plots; if x and or y is an
  ExpressionSet, then the first feature will be extracted
  from each; use the '[' operator to subset either x or y
  if necessary.
}
\examples{
\dontrun{
	hent3 <- x["SLC29A3", ]
	dotplot(hent3, sort=T, samples="APGI_1966")

	dotplot(x, feature="SLC29A3", sort=T, samples="APGI_1966")

 num <- rnorm(80); names(num) <- sampleNames(x)[1:80]
 dotplot(hent3, num)
 dotplot(num, hent3)

 dotplot(x["BRCA2",], x["BRCA1", ])
 dotplot(x, x, feature="BRCA1")
}
}
\author{
  Mark Cowley, 2012-10-16
}

