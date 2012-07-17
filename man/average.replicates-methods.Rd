\docType{methods}
\name{average.replicates}
\alias{average.replicates}
\alias{average.replicates,AnnotatedDataFrame,character-method}
\alias{average.replicates,AnnotatedDataFrame,factor-method}
\alias{average.replicates,ANY,missing-method}
\alias{average.replicates,data.frame,ANY-method}
\alias{average.replicates,eSet,character-method}
\alias{average.replicates,eSet,factor-method}
\alias{average.replicates,LumiBatch,character-method}
\alias{average.replicates,LumiBatch,factor-method}
\alias{average.replicates,matrix,character-method}
\alias{average.replicates,matrix,factor-method}
\title{Average replicate columns}
\arguments{
  \item{x}{a \code{matrix}, \code{data.frame}, \code{eSet},
  \code{LumiBatch} or \code{AnnotatedDataFrame} of data}

  \item{classes}{a \code{character} or \code{factor}: the
  sample classes, 1 per sample in \code{x}, with as many
  unique values or \code{levels} as there are unique
  samples.}
}
\value{
  The same data type as \code{x}, but with fewer samples,
  where the replicates have been averaged.
}
\description{
  Average the replicate samples within a dataset \code{x},
  based on the groupings indicated in by \code{classes},
  where identical values indicate samples to be
  \code{averaged}. Note that \code{averaged} is a bit of a
  generalisation, since some dataslots, eg \code{detection}
  within \code{LumiBatch} objects contains p-values, which
  are not averaged per se, see the
  \code{average.pval.replicates} S4 methods.
}
\details{
  \code{x} can be either a \code{matrix},
  \code{data.frame}, \code{eSet}, \code{LumiBatch},
  \code{AnnotatedDataFrame}. Note that for most of these
  data types, each sample is assumed to be in a column.
  \code{AnnotatedDataFrame} has 1 row per sample, thus the
  result will have fewer rows than \code{x}.

  \code{classes} can be a character or factor vector. Note
  that the \code{sampleNames} of the result are derived
  from the unique values, or levels in the \code{classes}
  vector

  \code{matrix} will return a \code{matrix}, with 1 column
  per \code{levels(classes)}, in the same order as
  \code{classes}. Note this is different

  \code{eSet} is a virtual object, thus all sub-classes of
  eSet, like ExpressionSet or AffyBatch should end up
  calling this method. This method will work as long as
  there are no additional slots defined in the class (which
  is why I have a separate LumiBatch-specific method).
  \code{average.replicates} on \code{eSet}'s does: (1)
  resizes \emph{all} elements within \code{AssayData},
  (which should all be matrix-like, and includes exprs,
  se.exprs & any other items found), as well as (2)
  \code{phenoData} and (3) \code{protocolData}, which
  indeed are the only elements defined in the \code{eSet}
  interface, which depend on the number of samples.\cr IF
  you get an error indicating that not all slots in the
  result contain the correct number of samples, then you
  may need to create another S4 method for that particular
  \code{eSet-sub-class}. Use the method implementation for
  \code{eSet} as a starting point, and then the method for
  \code{LumiBatch} as a worked example of an
  \code{eSet-sub-class} with a few additional slots. For
  (2,3) above, which are AnnotatedDataFrame's the replicate
  samples are averaged using the \code{average.replicates}
  for AnnotatedDataFrame method.

  \code{\linkS4class{LumiBatch}} objects contain a few
  additional slots in addition to those required by the
  \code{\linkS4class{eSet}} interface: \code{QC},
  \code{controlData}. Also, one of the \code{assayData}
  slots, \code{detection} contains detection p-values which
  shouldn't really just be averaged. This uses
  \code{average.pval.replicates} to properly average the
  detection pvalues from \code{\linkS4class{LumiBatch}}
  objects with replicates.

  \code{AnnotatedDataFrame} contains 1 row per sample,
  there's nothing in the spec as to required columns. The
  default for \code{LumiBatch} objects is a single column
  called sampleID, but this does not have to be the case.
  It's not clear how to average \code{character} values, so
  we retain the rows for the first sample matching each
  unique \code{class}.
}
\examples{
\dontrun{
load("Rmisc/x.norm.RDa.gz")
classes <- sub("\\\\.[12]$", "", sampleNames(x.norm))

# Matrix example
mat <- exprs(x.norm)
avg <- average.replicates(mat, classes)

# data.frame example
df <- as.data.frame(mat)
df.avg <- average.replicates(df, classes)

# LumiBatch example
x.av <- average.replicates(x.norm, classes)

# AnnotatedDataFrame example
adf <- PhenoData(x.norm)
adf.avg <- average.replicates(adf, classes)
adf.avg
as(adf.avg, "data.frame")
}
}
\author{
  Mark Cowley, 2011-09-01
}
\seealso{
  \code{\link{average.pval.replicates}}
}

