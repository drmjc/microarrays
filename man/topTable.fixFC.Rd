\name{topTable.fixFC}
\alias{topTable.fixFC}
\title{'fix' the logFC column for biologists (sorry!)}
\usage{
  topTable.fixFC(tt, digits = 4)
}
\arguments{
  \item{tt}{the output from \code{\link[limma]{topTable}}}

  \item{digits}{how many digits? see \code{\link{round}}}
}
\value{
  an edited topTable, with 3 extra columns (eg for
  logFC=-1.82):\cr \dQuote{FC}: the unlogged fold change
  (eg 0.283)\cr \dQuote{direction}: up or down (eg down)\cr
  \dQuote{absFC}: the unlogged absolute fold change (eg
  3.531)\cr
}
\description{
  'fix' the logFC column for biologists (sorry!)
}
\examples{
\dontrun{
head( topTable.fixFC(tt) )
#            ID logFC     t  P.Value adj.P.Val     B    FC direction absFC
# 16143 7923547  2.58 13.44 3.31e-06    0.0695 -1.41 5.985        up  5.99
# 12361 8096301  2.59  9.52 3.22e-05    0.3382 -1.62 6.039        up  6.04
# 16282 7953200  1.74  8.62 6.11e-05    0.4273 -1.70 3.331        up  3.33
# 540   7995729  1.64  8.03 9.52e-05    0.4996 -1.77 3.109        up  3.11
# 8641  8141374 -1.42 -7.43 1.55e-04    0.6164 -1.85 0.374      down  2.67
# 11546 8007848  1.48  7.28 1.76e-04    0.6164 -1.87 2.790        up  2.79
}
}
\author{
  Mark Cowley, 2008-07-14
}

