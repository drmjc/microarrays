#' @name as.data.frame.AnnotatedDataFrame
#' @title Coerce to a Data Frame
#' @method as.data.frame AnnotatedDataFrame
#' @importClassesFrom Biobase AnnotatedDataFrame
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
