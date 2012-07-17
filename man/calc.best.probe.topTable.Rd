\name{calc.best.probe.topTable}
\alias{calc.best.probe.topTable}
\title{From a toptable which contains probe ID's in the ID column, determine the
best probe for each gene}
\usage{
  calc.best.probe.topTable(tt, probe2gene, toupper = TRUE,
    symbols.ignore.list = c("---", ""), only.genes = FALSE)
}
\arguments{
  \item{tt}{a toptable}

  \item{probe2gene}{a 2 column data frame with probes and
  gene symbols}

  \item{toupper}{if TRUE, convert the gene symbol to
  UPPERCASE; FALSE leaves it untouched}

  \item{symbols.ignore.list}{Probes with these values as
  gene symbols, as well as NA and NULL will not be merged
  to the same value, otherwise the 1000 probes that map to
  "---" will be collapsed to a single 'gene'}

  \item{only.genes}{After collapsing multiple probes to
  one, if TRUE, then only the probes with real gene symbols
  will be exported; if FALSE, all collapsed probes with or
  without a gene symbol will be exported.}
}
\description{
  If multiple probes map to a gene symbol, then just 1 row
  will be returned - the one with the max(abs(tstat)) If
  the probe has no symbol (or is "---", "", NA, NULL,
  symbols.ignore.list) then it will also be exported.
}
\author{
  Mark Cowley, 2008-12-08
}

