\name{merge_topTable}
\alias{merge_topTable}
\title{Function to merge N topTable's together, retaining a particular column from
each topTable.}
\usage{
  merge_topTable(tt.list,
    keep = c("t", "adj.P.Val", "P.Value"))
}
\arguments{
  \item{tt.list}{a named list of topTable objects.}

  \item{keep}{the column to keep or retain. It should be
  one of the column names in each GSEA object.}
}
\value{
  a data.frame of N+1 columns: ID, then the values that
  were extracted from each topTable. eg the t-stats
}
\description{
  Function to merge N topTable's together, retaining a
  particular column from each topTable.
}
\author{
  Mark Cowley, 2009-12-16
}

