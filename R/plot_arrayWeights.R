#' Plot array weights
#' 
#' Function to plot (and optionally calculate) the array weights using
#' \code{\link[limma]{arrayWeights}} in \code{limma}.
#' 
#' @param aw a named, numeric vector of array weights [optional]
#' @param data if \code{aw=NULL}, then you can specify an object that must be
#'   compatible with arrayWeights. see \code{\link[limma]{arrayWeights}}
#' @param design if \code{aw=NULL}, you can specify the design matrix for calculating
#'   the array weights. defaults to the unit vector.
#' @param verbose logical
#' @param main the plot title
#' @param \dots arguments passed to \code{\link{barplot}}
#' @return creates a barplot of array weights, and invisibly returns the array weights.
#' @author Mark Cowley, 2008-07-11
#' @export
#' @importFrom limma arrayWeights
plot_arrayWeights <- function(aw=NULL, data=NULL, design=NULL, verbose=TRUE, main="Array Weights", ...) {
	if( is.null(data) && is.matrix.like(aw) ) {
		stop("Have you forgotten to specfic the data= argument")
	}
	if( is.null(aw) ) {
		if( is.null(design) )
			design <- rep(1, ncol(data))
		if( verbose )
			cat("Fitting arrayWeights linear model.\n")
		aw <- arrayWeights(as.matrix(data), design=design)
		names(aw) <- colnames(data)
	}
	opar <- par(no.readonly=TRUE)
	par(las=2)
	barplot(aw, names=names(aw), main=main, ylab="weight", ...)
	abline(h=1, lty="dashed")
	par(opar)
	
	invisible( aw )
}
