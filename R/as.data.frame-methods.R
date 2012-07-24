#' @name as.data.frame.AnnotatedDataFrame
#' @title Coerce to a Data Frame
#' @S3method as.data.frame AnnotatedDataFrame
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @aliases as.data.frame
#' @usage as.data.frame(x, ...)
#' 
#' @param x an AnnotatedDataFrame
#' @param \dots unused
#' 
#' @return a \code{data.frame}
#' 
#' @author Mark Cowley, 2012-07-24
as.data.frame.AnnotatedDataFrame <- function(x, ...) {
	as(x, "data.frame")
}
