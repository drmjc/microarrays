\name{plot_cor_topTable}
\alias{plot_cor_topTable}
\title{Compare two topTables via correlation of the t-stats}
\usage{
  plot_cor_topTable(tt1, tt2, xlab = "topTable 1",
    ylab = "topTable 2", main = "topTable correlation",
    column = "t", ...)
}
\arguments{
  \item{tt1}{two instances of topTables. They can have
  different numbers of probes, but some overlap is
  required.}

  \item{tt2}{two instances of topTables. They can have
  different numbers of probes, but some overlap is
  required.}

  \item{xlab}{obvious}

  \item{ylab}{obvious}

  \item{main}{obvious}

  \item{column}{which column to plot. default is "t".
  "logFC" would be another good choice.}

  \item{\dots}{additional arguments passed to plot.cor. see
  plot.cor}
}
\value{
  creates a correlation plot. see plot.cor invisibly
  returns the 2 toptables merged on a common ID. The unique
  suffixes on the colnames comes from xlab & ylab.  See
  Also: plot.cor
}
\description{
  Compare two topTables via correlation of the t-stats
}
\author{
  Mark Cowley, 2010-01-18
}

