\name{plot_genespringGX_lines}
\alias{plot_genespringGX_lines}
\title{Perform a GeneSpring GX 7.3.1 line plot}
\usage{
  plot_genespringGX_lines(data,
    row.centre = c("median", "mean", "none"),
    colour.by = 1, symmetrical = TRUE, main = "",
    xlab = "", ylab = "Expression Ratio (log2)", las = 2,
    colour.by.line = TRUE, lwd = 0.5, bg.col = "black",
    alpha = 0.2, ...)
}
\arguments{
  \item{data}{\code{matrix} or \code{data.frame} of gene
  expression data}

  \item{row.centre}{the method to use to centre each gene
  expression trait. one of \dQuote{median}, \dQuote{mean},
  or \dQuote{none}}

  \item{colour.by}{which array to use to set the colour}

  \item{symmetrical}{logical: make the y-axis symmetrical}

  \item{main}{the plot title}

  \item{xlab}{the x-axis title}

  \item{ylab}{the y-axis title}

  \item{las}{see \code{\link{par}}. default=2}

  \item{colour.by.line}{logical: if \code{TRUE} add a
  vertical line to indicate which array was used to set the
  colour (see \code{colour.by})}

  \item{lwd}{line width}

  \item{bg.col}{the default background colour.
  default=\dQuote{black}}

  \item{alpha}{control the opacity of the lines.
  default=0.2}

  \item{\dots}{arguments passed to
  \code{\link{plot.matrix}}}
}
\value{
  none.
}
\description{
  Perform a GeneSpring GX 7.3.1 line plot
}
\examples{
m <- matrix(rnorm(1000*10),1000,10)
plot_genespringGX_lines(m)
}
\author{
  Mark Cowley
}

