#' subset an eSet object
#' subset the rows (features) and columns(samples) within an eSet object.
#' @note if you have a subclass of eSet, there may be slots that are not
#' subsetted. For example LumiBatch objets can have an extra controlData
#' slot.
#' @method subset eSet
#' @param x an eSet
#' @param subset logical expression indicating elements or rows to keep:
#'          missing values are taken as \code{FALSE}.
#' @param select logical expression, indicating columns to select from a \code{data.frame}.
#' @param drop logical: ignored
#' @param \dots ignored
#' @return an eSet with fewer rows and or columns
#' @author Mark Cowley, 2011-09-01
#' @seealso \code{\link{subset.LumiBatch}}
#' @export
#' @S3method subset eSet
subset.eSet <- function(x, subset, select, drop = FALSE, ...) {
	if( missing(subset) ) subset <- rep(TRUE, nrow(x))
	if( missing(select) ) select <- rep(TRUE, ncol(x))

	# The [ and ] operators work well for eSet's already.
	x[which(subset), which(select)]
}
# CHANGELOG
# 2011-09-01: v1
# 2011-09-22: bug fix, where subset was used for columns, and select for rows, instead of the opposite
# in the standard subset function.
# 



#' subset a LumiBatch object
#' 
#' This is preferable to subset.eSet, since there are extra slots in a 
#' LumiBatch
#' 
#' The support for subsetting the controlData is still inadequately tested.
#' For instance, do all controlData slots have the first 2 columns being
#' "?" and "ProbeID"?
#' 
#' @method subset LumiBatch
#' @param x an LumiBatch
#' @param subset logical expression indicating elements or rows to keep:
#'          missing values are taken as \code{FALSE}.
#' @param select logical expression, indicating columns to select from a \code{data.frame}.
#' @param drop logical: ignored
#' @param \dots ignored
#' @return an LumiBatch with fewer rows and or columns
#' @author Mark Cowley, 2011-09-01
#' @export
#' @S3method subset LumiBatch
subset.LumiBatch <- function(x, subset, select, drop=FALSE, ...) {
	if( missing(subset) ) subset <- rep(TRUE, nrow(x))
	if( missing(select) ) select <- rep(TRUE, ncol(x))

	# The [ and ] operators work well for eSet's already.
	res <- x[which(subset), which(select)]
	res@controlData <- subset(x@controlData,select=c(TRUE, TRUE, select))
	res
}
