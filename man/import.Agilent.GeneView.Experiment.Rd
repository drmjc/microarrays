\name{import.Agilent.GeneView.Experiment}
\alias{import.Agilent.GeneView.Experiment}
\title{Import a collection of GeneView files}
\usage{
  import.Agilent.GeneView.Experiment(files, names,
    debug = FALSE)
}
\arguments{
  \item{files}{the file names. (These will NOT be
  reordered)}

  \item{names}{if specified, then these names will be used
  in the column names of the data.frames; if NULL, then a
  name based on the filename will be used}

  \item{debug}{logical}
}
\value{
  a list with 5 elements as indicated in Details.
}
\description{
  Import a set of Agilent GeneView files that logically
  form an experiment.
}
\details{
  This first imports each GeneView file, the collates the
  information into 3 \code{data.frame}s, with ngenes x
  nsamples, called:\cr \code{value}: the gTotalGeneSignal
  column\cr \code{error}: the gTotalGeneError column\cr
  \code{detected}: the gIsDetected column\cr In addition to
  these, there is a \code{vector} from the
  \code{ControlType} column (which is identical across all
  arrays), and a \code{vector} called \code{ProbeIDs}.
}
\author{
  Mark Cowley, 2008-08-01
}

