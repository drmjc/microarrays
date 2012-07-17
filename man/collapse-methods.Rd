\docType{methods}
\name{collapse}
\alias{collapse}
\alias{collapse,LumiBatch,logical,logical,function,character-method}
\alias{collapse,LumiBatch,missing,missing,missing,missing-method}
\title{collapse}
\arguments{
  \item{x}{an object}

  \item{decreasing}{logical: should the sort order be
  increasing or decreasing?}

  \item{na.last}{for controlling the treatment of
  \code{NA}s.  If \code{TRUE}, missing values in the data
  are put last; if \code{FALSE}, they are put first; if
  \code{NA}, they are removed.}

  \item{FUN}{a function used to determine which row to
  include, when there are multiple rows with the same
  'column' value. eg \code{mean, min, max, var, sd, mad,
  prod}}

  \item{column}{the column name within \code{x} which
  contains the keys. see details}

  \item{\dots}{arguments passed to FUN}
}
\value{
  something
}
\description{
  collapse a dataset to have just 1 row per unique value of
  key.
}
\details{
  given an object with data, and 1 'key' column, containing
  possibly non-unique identifiers, create a result, which
  has unique values for those keys. The 'key' column could
  be a gene symbol, and each row could be a probe; the goal
  is to convert a 1-row-per-probe table into
  1-row-per-gene.
}
\section{collapse,LumiBatch and column}{
  The \code{column} should be the name of a column found in
  the \code{fData(x)} slot. hint: \code{colnames(fData(x)),
  or fvarLabels(x)}
}

\section{collapse,LumiBatch and Missing values}{
  Often some probes don't have a genesymbol, and should
  thus have an \code{NA} in this column. \code{na.list}
  controls what to do with the gene-less probes.
  \code{na.last=TRUE} keeps the probes & moves them to the
  bottom of the result. \code{na.last=FALSE} keeps the
  probes & leaves them in the sorted order specified by
  \code{FUN}. Note this is different to other
  implementations of \code{na.last}. \code{na.last=NA
  (default)} discards these probes from the result.
}
\examples{
\dontrun{
 load("Rmisc/x.averaged.RDa.gz")
	tmp <- collapse(x.averaged, FUN=var, decreasing=TRUE, na.last=FALSE, "SymbolReannotated")
}

require(lumi)
data(example.lumi)
featureData(example.lumi)$GeneSymbol <- c("TP53", "BRCA1", "BRCA2", "KRAS", "SMAD4")
collapse(example.lumi, TRUE, NA, mean, "GeneSymbol")
}
\author{
  Mark Cowley
}

