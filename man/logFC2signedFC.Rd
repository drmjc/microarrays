\name{logFC2signedFC}
\alias{logFC2signedFC}
\title{Convert logFC (FoldChange) to signed FC
When logFC<0, unlogged FC values get squashed in (0,1). The signed FC
has the same magnitude whether the logFC was >0 or <0, but if <0, a leading - is
added}
\usage{
  logFC2signedFC(fc)
}
\arguments{
  \item{fc}{a numeric vector of log-base-2 Fold Change
  values}
}
\value{
  a numeric vector of signed fold change values
}
\description{
  Convert logFC (FoldChange) to signed FC When logFC<0,
  unlogged FC values get squashed in (0,1). The signed FC
  has the same magnitude whether the logFC was >0 or <0,
  but if <0, a leading - is added
}
\author{
  Mark Cowley, 2009-02-04
}

