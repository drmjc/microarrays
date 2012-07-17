\name{uncollapse.rows}
\alias{uncollapse.rows}
\title{Opposite of collapse.rows.}
\usage{
  uncollapse.rows(x, cols2uncollapse = NULL, sep = " // ")
}
\arguments{
  \item{cols2uncollapse}{Which columns need uncollapsing?
  Must specify at least 1 column (hint: this is the column
  that contains \code{sep} that you're trying to get rid
  of). If you specify >1 columns, then each cell in that
  row must have the same number of code elements to be
  split.}

  \item{x}{a \code{data.frame}}

  \item{sep}{the seperator, such as \dQuote{, } or \dQuote{
  // }}
}
\description{
  Strongly suggest using this function to reverse the
  effects of \code{\link{collapse.rows}}, using the same
  arguments that were supplied to
  \code{\link{collapse.rows}} itself.
}
\details{
  Collapsed data means a \code{data.frame} with at least 1
  column whose values are \code{sep} delimited. Eg a row of
  data with the gene symbol "Ankrd11|Galnt2" or a row of
  data with multiple GO terms, eg "GO:00001 /// GO:00002
  /// GO:00003". This function takes this data, and
  increases the number of rows, such that these data have 1
  element per row. so "Ankrd11" and "Galnt2" for example.
  Thus changing it from n:1 to 1:1.

  All columns that are not specified in
  \code{cols2uncollapse} will be repeated. If you have just
  1 column to uncollapse, then only that column will be
  changed. If you have more than 1 column to expand, then
  within those rows that need uncollapsing, all specified
  columns MUST have the same number of elements. Eg
  consider a \code{data.frame} with 1 row per gene with 3
  GO-term columns: GO.ID, GO.Name, GO.Evidence. For any
  given gene with 3 GO terms, there should also be exactly
  # GO ID's, 3 GO Names and 3 GO term evidence codes. If
  there are different numbers of elements found this will
  thow an error.
}
\examples{
df <- data.frame(
   Name=rep(LETTERS[1:3], each=3),
   Description=rep(letters[1:3], each=3),
   Value=LETTERS[11:19],
   stringsAsFactors=FALSE
)
a <- collapse.rows(df, 1, 3)
uncollapse.rows(a, 1, 3)
}
\author{
  Mark Cowley, 2009-01-08
}
\seealso{
  \code{\link{collapse.rows}}
}

