#' dotplots on ExpressionSet objects
#' 
#' Plot either a single feature from an ExpressionSet as a dotplot (default is sorted low to high), or
#' generate a 2D dotplot from combinations of 2 ExpressionSet's, or ExpressionSet and a named numeric, where
#' those names overlap with the featureNames from the ExpressionSet. As always, the 2D objects need to
#' have some overlap in their names, but otherwise, data unique to each set will be
#' silently ignored.
#' 
#' @section 2D dotplots:
#' these are just xy plots; if x and or y is an ExpressionSet, then the first feature will
#' be extracted from each; use the '[' operator to subset either x or y if necessary.
#'
#' @param x an ExpressionSet
#' @param y an ExpressionSet, or missing.
#' @param sort logical: if \code{TRUE}, then sort from low to high (default=TRUE). only applies if y is missing.
#' @param feature the index or name of the feature to plot. default=1. In the 2D case, if the feature exists in
#' both ExpressionSet's then that will be seleected. If plotting different features from x and y, then leave feature=1,
#' and select the feature of interest from each expression set using '['. see examples.
#' @param xlab see par
#' @param ylab see par
#' @param main see par
#' @param add.mean logical: add a horizontal line about the mean of that feature
#' @param samples an optional character vector of sample names to highlight and label in red.
#' @param lowess.col the colour of the 2D loess line. \code{NA} to not show one, otherwise choose a named colour.
#' @return nothing
#' 
#' @author Mark Cowley, 2012-10-16
#' 
#' @importFrom Biobase featureNames sampleNames
#' @importMethodsFrom mjcgraphics dotplot
#' @export
#' 
#' @rdname dotplot-methods
#' @aliases dotplot,ExpressionSet,missing,ANY-method
#' 
#' @examples
#' \dontrun{
#' 	hent3 <- x["SLC29A3", ]
#'	dotplot(hent3, sort=T, samples="APGI_1966")
#' 
#'	dotplot(x, feature="SLC29A3", sort=T, samples="APGI_1966")
#' 
#'  num <- rnorm(80); names(num) <- sampleNames(x)[1:80]
#'  dotplot(hent3, num)
#'  dotplot(num, hent3)
#' 
#'  dotplot(x["BRCA2",], x["BRCA1", ])
#'  dotplot(x, x, feature="BRCA1")
#' }
setMethod("dotplot",
	signature=signature("ExpressionSet", "missing", "ANY"),
	function(x, y, sort, feature=1, xlab="Rank", ylab="Expression Level (log2)", main=featureNames(x)[feature], add.mean=TRUE, samples=NULL, ...) {
		if( missing(sort) ) sort <- TRUE
		o <- 1:ncol(x)
		if( sort ) {	
			o <- order(exprs(x)[feature,,drop=T], decreasing=FALSE)
		}
		dotplot(exprs(x)[feature,o,drop=T], xlab="Rank", ylab="Expression Level (log2)", main=main, sort=sort, ...)
		if( add.mean ) {
			abline(h=rowMeans(exprs(x[feature, ])), col="red")
			text(1, rowMeans(exprs(x[feature, ])), "mean", col="red", pos=3)
		}
		if( !is.null(samples) ) {
			if( !all(samples %in% sampleNames(x)) ) stop("some samples not in sampleNames(x): ", setdiff(samples, sampleNames(x)))
			idx <- which(sampleNames(x) == samples)
			ypos <- exprs(x[feature,idx])[1,,drop=T]
			xpos <- which(o == idx)
			points(xpos, ypos, pch=19, col="red")
			text(xpos, ypos, samples, col="red", pos=2)
		}
	}
)

# dotplot.ExpressionSet <- function(x, feature=1, xlab="Rank", ylab="Expression Level (log2)", main=featureNames(x)[feature], add.mean=TRUE, sort=TRUE, samples=NULL) {
# 	o <- 1:ncol(x)
# 	if( sort ) {	
# 		o <- order(exprs(x)[feature,,drop=T], decreasing=FALSE)
# 	}
# 	plot(exprs(x)[feature,o,drop=T], xlab="Rank", ylab="Expression Level (log2)", main=main)
# 	if( add.mean ) {
# 		abline(h=rowMeans(exprs(x[feature, ])), col="red")
# 		text(1, rowMeans(exprs(x[feature, ])), "mean", col="red", pos=3)
# 	}
# 	if( !is.null(samples) ) {
# 		if( !all(samples %in% sampleNames(x)) ) stop("some samples not in sampleNames(x): ", setdiff(samples, sampleNames(x)))
# 		xpos <- o[which(sampleNames(x) == samples)]
# 		ypos <- exprs(x[feature,xpos])[1,,drop=T]
# 		points(xpos, ypos, pch=19, col="red")
# 		text(xpos, ypos, samples, col="red", pos=2)
# 	}
# }


#' @rdname dotplot-methods
#' @aliases dotplot,ExpressionSet,ExpressionSet,ANY-method
#' @importFrom stats lowess
setMethod("dotplot",
	signature=signature("ExpressionSet", "ExpressionSet", "ANY"),
	function(x, y, sort, feature=1, loess.col=NA, ...) {
		common.names <- intersect(sampleNames(x), sampleNames(y))
		length(common.names) > 0 || stop("no sample names found in common")
		x <- x[,common.names]
		y <- y[,common.names]
		a <- exprs(x[feature, ])[1,,drop=T]
		b <- exprs(y[feature, ])[1,,drop=T]
		names(a) <- names(b) <- common.names
		
		has.data <- !is.na(a) & !is.na(b)
		
		dotplot(a, b, ...)
		title(sub=paste("n =", sum(has.data)), outer=FALSE)
		
		if( !is.na(loess.col) ) {
			lines(stats::lowess(a[has.data],b[has.data]), col=loess.col)
		}
	}
)


#' @rdname dotplot-methods
#' @aliases dotplot,ExpressionSet,numeric,ANY-method
setMethod("dotplot",
	signature=signature("ExpressionSet", "numeric", "ANY"),
	function(x, y, sort, feature=1, ...) {
		x <- exprs(x[feature,])[1,,drop=T]
		dotplot(x, y, ...)
	}
)


#' @rdname dotplot-methods
#' @aliases dotplot,numeric,ExpressionSet,ANY-method
setMethod("dotplot",
	signature=signature("numeric", "ExpressionSet", "ANY"),
	function(x, y, sort, feature=1, ...) {
		y <- exprs(y[feature,])[1,,drop=T]
		dotplot(x, y, ...)
	}
)
