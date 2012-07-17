\name{export.topTable}
\alias{export.topTable}
\title{Export a topTable, or a list of top tables to an excel file.
Export a topTable, or a list of top tables to an excel file, with an
optional summary of the number of DE genes passing various generic
thresholds.}
\usage{
  export.topTable(tt, file, annot = NULL, fixFC = TRUE,
    Pcount = NULL, coefficients = NULL, fit = NULL,
    digits = 4, summary = TRUE, drop.Bstat = TRUE,
    adj.Pval.colname = "FDR", ...)
}
\arguments{
  \item{tt}{a topTable data.frame}

  \item{file}{the output file name}

  \item{annot}{the probe-level annotation object. rownames
  should be set to the probe ID, so tt$ID can match them}

  \item{fixFC}{fix the logFC to produce an absFC and
  direction columns}

  \item{Pcount}{an optional vector of Pcounts (ie how many
  times each probe was detected as present). Should be
  named so that tt$ID can match the probes}

  \item{coefficients}{an optional data.frame, of
  coefficients from the lmFit.  You can use this arg to add
  custom columns to the topTable. eg the AvgExpr}

  \item{fit}{the lmFit object. Used to obtain the Standard
  error of each probe. using fit$stdev.unscaled *
  fit$sigma}

  \item{digits}{how many digits to round to?}

  \item{summary}{insert a worksheet into the excel workbook
  which contains a summary of the toptable}

  \item{drop.Bstat}{logical. If TRUE then the "B" column is
  not exported.}

  \item{adj.Pval.colname}{rather than "adj.P.Val", you can
  rename to column to "FDR", or "q", or similar.}

  \item{\dots}{additional arguments passed to write.xls}
}
\description{
  The first worksheet will be the summary, and subsequent
  worksheets are for each top table supplied. Additional
  benefits:\cr - Limits the number of significant figures
  to digits correctly for the P and q values\cr - includes
  a Pcount column\cr - includes gene information which uses
  the ID column and rownames of annot object\cr - includes
  a few columns of coefficients, for eg if the group-means
  parameterisation is used, then the values of the
  coefficients (ie the value for each 'group') may be
  appropriate to export.
}
\author{
  Mark Cowley
}

