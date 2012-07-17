\name{topTable.SE}
\alias{topTable.SE}
\title{Add a moderated SE column to a topTable object.}
\usage{
  topTable.SE(fit, coef = NULL, ...)
}
\arguments{
  \item{fit}{see topTable}

  \item{coef}{the coefficient index. if \code{NULL}, this
  is set to 1.}

  \item{\dots}{arguments passed to
  \code{\link[limma]{topTable}}}
}
\value{
  a topTable object produced by
  \code{\link[limma]{topTable}}, in addition to an
  \dQuote{SE} column containing the moderated standard
  errors.
}
\description{
  moderated SE's are fit$stdev.unscaled * sqrt(fit$s2.post)
  unmoderated SE's are fit$stdev.unscaled * fit$sigma from
  p53 & 54 limma userguide: for gene j, contrast k: ujk =
  fit$stdev.unscaled = unscaled stdev sj = fit$sigma = gene
  j's residual variance, with dj df Bjk = linear model
  estimates = the logFC's tstat.ord <- signal/SE <-
  Bjk/(ujk.sj) <- fit$coef/(fit$stdev.unscaled * fit$sigma)
  tstat.mod <- signal/SE' <- Bjk/(ujk.s'j) <-
  fit$coef/(fit$stdev.unscaled * sqrt(fit$s2.post)) ^^^^
  which follows a t-dist with d0 + dj d.f. s2.post is the
  weighted avg of s2.prior and sigma^2 with weights
  proportional to df.prior and df.residual, respectively.
}
\author{
  Mark Cowley, 2010-11-12
}

