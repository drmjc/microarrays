#' Perform a GeneSpring GX 7.3.1 line plot
#' @param data \code{matrix} or \code{data.frame} of gene expression data
#' @param row.centre the method to use to centre each gene expression trait. one of \dQuote{median},
#' \dQuote{mean}, or \dQuote{none}
#' @param colour.by which array to use to set the colour
#' @param symmetrical logical: make the y-axis symmetrical
#' @param main the plot title
#' @param xlab the x-axis title
#' @param ylab the y-axis title
#' @param las see \code{\link{par}}. default=2
#' @param colour.by.line logical: if \code{TRUE} add a vertical line to indicate which array
#'   was used to set the colour (see \code{colour.by})
#' @param lwd line width
#' @param bg.col the default background colour. default=\dQuote{black}
#' @param alpha control the opacity of the lines. default=0.2
#' @param \dots arguments passed to \code{\link{plot.matrix}}
#' @return none.
#' @author Mark Cowley
#' @export
#' @examples
#' m <- matrix(rnorm(1000*10),1000,10)
#' plot_genespringGX_lines(m)
plot_genespringGX_lines <- function(data, row.centre=c("median", "mean", "none"), colour.by=1, symmetrical=TRUE,
	main="", xlab="", ylab="Expression Ratio (log2)", las=2, colour.by.line=TRUE, lwd=0.5, bg.col="black", alpha=0.2, ...) {
	opar <- par(no.readonly=TRUE)
	on.exit(par(opar))
	par(las=las)
	
	row.centre <- row.centre[1]
	if( row.centre == "median" )
		vals <- apply(data, 1, median)
	else if( row.centre == "mean" )
		vals <- apply(data, 1, mean)
	else if( row.centre == "none" )
		vals <- rep(0, nrow(data))
	data <- data - vals

	cols <- colour.step(from="blue", to="red", via="yellow", steps=nrow(data), alpha=alpha)[rank(data[, colour.by])]

	ylim <- range(data, na.rm=TRUE)
	if( symmetrical )
		ylim <- symmetricise(ylim)

	plot.matrix(data, main=main, xlab=xlab, ylab=ylab, col=cols, xlabels=colnames(data), ylim=ylim, lwd=lwd, bg.col=bg.col, ...)
	if( colour.by.line ) abline(v=colour.by, lty="dotted", col="white")
}
