\name{collapse.topTable}
\alias{collapse.topTable}
\title{Take a topTable with one row per probeset, and collapse it to one row per
gene symbol.}
\usage{
  collapse.topTable(tt, probe2gene, toupper = FALSE,
    symbols.ignore.list = c("---", ""), only.genes = FALSE,
    verbose = TRUE)
}
\arguments{
  \item{tt}{a toptable, or a list of toptables}

  \item{probe2gene}{a 2 column data.frame, 1 row per probe,
  mapping probe ID's (column 1) to gene symbols (column 2)}

  \item{toupper}{convert the gene symbol in column 2 to
  UPPER case?}

  \item{symbols.ignore.list}{a vector of symbols that are
  used to identify probes with no gene. NB: the special
  values NA and NULL are always used, so no need to specify
  these.}

  \item{only.genes}{if TRUE, then filter out the probes
  without a gene symbol.  if FALSE, then probes with or
  without gene symbols are exported.}

  \item{verbose}{Helpful messages?}
}
\value{
  a data.frame like tt, but with some new columns:
  "Gene.Symbol" in pos 2, "Nprobes" = the number of
  possible probes for the given gene, "Possible.Probes" = a
  comma separated character(1) of the possible probes for
  the given gene. nrows = number of unique genes +
  if(only.genes=FALSE) the number of probes with no gene.
  row ordering is given by the original tt's row ordering.
}
\description{
  Probes either map to Gene Symbols, or not & some Genes
  can have multiple probes mapping to them (due to array
  design, or isoforms). This method will choose the best
  probe per gene, by selecting that with the largest
  abs(t.stat), since that probe 'performed' the best, (can
  be either up or down-regulated). In my experience, for
  those genes with lots of probes, many do not hyb well (eg
  on the HGU133+2 array, due to poor probe design), and one
  will hyb better than the others & thus produce larger +/-
  t.stat. This is the best peforming probe & should be
  selected to represent that gene. Other alternatives
  include median probe, max.avg, ... (not yet implemented.)
  Probes that don't map to symbols can either be filtered
  out (only.genes=TRUE), or left in (only.genes=FALSE).
  Probes that don't map to gene symbols are identified via
  the probe2gene table & either have NA, NULL, or
  symbols.ignore.list in the 2nd column of the probe2gene.
  You can add custom values to the symbols.ignore.list, eg
  c("N/A", "missing", "blank", "neg.con", ...)
}
\author{
  Mark Cowley, 2011-02-22
}

