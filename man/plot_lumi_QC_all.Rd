\name{plot_lumi_QC_all}
\alias{plot_lumi_QC_all}
\title{Perform all of the QC plots that are mentioned in the lumi vignette.
It creates a density plot, boxplot, PCA, HCL, CV-plot, pairs and MA plots.
The pairs and MA plot can take a LONG time (>10mins?) if you have lots of arrays,
so you can skip these plots using the \code{MA=FALSE} and \code{pairs=FALSE}
arguments}
\usage{
  plot_lumi_QC_all(x, dir, prefix,
    title.prefix = "Unnormalised", MA = TRUE, pairs = TRUE)
}
\description{
  Perform all of the QC plots that are mentioned in the
  lumi vignette. It creates a density plot, boxplot, PCA,
  HCL, CV-plot, pairs and MA plots. The pairs and MA plot
  can take a LONG time (>10mins?) if you have lots of
  arrays, so you can skip these plots using the
  \code{MA=FALSE} and \code{pairs=FALSE} arguments
}
\note{
  The CV plot some times fails, even if you don't have
  missing data in your LumiBatch. In this case, you'll get
  this Error:\cr \code{Error in density.default(newX[, i],
  ...) : 'x' contains missing values}
}
\examples{
\dontrun{
dir.create("QC/01.unnorm")
plot_lumi_QC_all(x.raw, "QC/01.unnorm/", "raw", "Unnormalised")
dir.create("QC/02.transformed")
plot_lumi_QC_all(x.transformed, "QC/02.transformed/", "vst", "VST
Transformed")
dir.create("QC/03.rsn")
plot_lumi_QC_all(x.norm.rsn, "QC/03.rsn/", "rsn", "RSN Normalised")
}
}
\author{
  Mark Cowley, 2008-10-23
}

