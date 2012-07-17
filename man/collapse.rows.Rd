\name{collapse.rows}
\alias{collapse.rows}
\title{collapse rows}
\usage{
  collapse.rows(x, key.column = 1, cols2collapse = NULL,
    sep = " // ")
}
\arguments{
  \item{x}{a \code{data.frame}}

  \item{key.column}{the column that must end up having one
  row per key. numeric or character allowed.}

  \item{cols2collapse}{which are the columns that you want
  to collapse. Often there are columns which will contain
  the same info repeated over and over, and you don't want
  these things to have the same word repeated N times.
  Just ignore these columns then, and only supply the
  column names of those columns that you want to be
  joined.}

  \item{sep}{the seperator, such as \dQuote{, } or \dQuote{
  // }}
}
\value{
  a \code{data.frame} with same num columns, but only N
  rows corresponding to the N different values in the key
  column. alphanumerically sorted by key column.
}
\description{
  Take a \code{data.frame} which has rows that contain
  mostly the same info, but some columns change. You want
  one row per 'key' which is in a particular column, and in
  the columns that contain non-equal data, collapse these
  values into a string, where each entry is separated by ",
  " or " // " for example.
}
\examples{
df <- data.frame(
   Name=rep(LETTERS[1:3], each=3),
   Description=rep(letters[1:3], each=3),
   Value=LETTERS[11:19],
   stringsAsFactors=FALSE
)
collapse.rows(df, 1, 3)
}
\author{
  Mark Cowley, 2009-01-07
}
\seealso{
  \code{\link{uncollapse.rows}}
}

