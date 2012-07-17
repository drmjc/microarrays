\name{plot_venn3D_topTable}
\alias{plot_venn3D_topTable}
\title{Venn Diagram from 3 topTable objects
Plot a Venn Diagram of the differentially expressed genes within 3 topTable objects.
Generates a 4x3 panel of plots, 1 per each p.threshold, q.thresh, logFC.thresh and
top N sizes, where N defaults to 50, 100, 200.}
\usage{
  plot_venn3D_topTable(tt1, tt2, tt3,
    p.thresh = c(0.05, 0.001, 1e-04),
    q.thresh = c(0.25, 0.1, 0.05),
    logFC.thresh = c(0.585, 1, 2), sizes = c(50, 100, 250),
    names = LETTERS[1:3])
}
\arguments{
  \item{tt1}{a \code{data.frame} from
  \code{\link[limma]{topTable}}}

  \item{tt2}{a \code{data.frame} from
  \code{\link[limma]{topTable}}}

  \item{tt3}{a \code{data.frame} from
  \code{\link[limma]{topTable}}}

  \item{p.thresh}{a vector of P Value thresholds, 1 per
  Venn Diagram}

  \item{q.thresh}{a vector of FDR thresholds, 1 per Venn
  Diagram}

  \item{logFC.thresh}{a vector of logFC thresholds, 1 per
  Venn Diagram}

  \item{sizes}{a vector of top N sizes}

  \item{names}{the names of the 3 topTable objects.
  default=\dQuote{A}, \dQuote{B}, \dQuote{C}}
}
\value{
  none.
}
\description{
  Venn Diagram from 3 topTable objects Plot a Venn Diagram
  of the differentially expressed genes within 3 topTable
  objects. Generates a 4x3 panel of plots, 1 per each
  p.threshold, q.thresh, logFC.thresh and top N sizes,
  where N defaults to 50, 100, 200.
}
\author{
  Mark Cowley, 2011-08-02
}

