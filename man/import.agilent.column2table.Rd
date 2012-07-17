\name{import.agilent.column2table}
\alias{import.agilent.column2table}
\title{Extract a column from Agilent files.
Take a vector of filenames, each of which is a large table from Agilent,
read each table in, then pull out a single column from each table, and
stitch these columns together into 1 large table.}
\usage{
  import.agilent.column2table(files, colname = "LogRatio")
}
\arguments{
  \item{files}{a vector of file names, including their
  paths}

  \item{colname}{the relevant column name to pull out.
  Beware of colnames with spaces.}
}
\value{
  a table of \code{N} columns, correspnding to the \code{N}
  files.
}
\description{
  Extract a column from Agilent files. Take a vector of
  filenames, each of which is a large table from Agilent,
  read each table in, then pull out a single column from
  each table, and stitch these columns together into 1
  large table.
}
\author{
  Mark Cowley, 2008-08-29
}

