#' sort
#'
#' @inheritParams base::sort
#' @param FUN a sort function. if \code{\link{x}} is 2D or more, then this is applied to the 1st
#' dimension (ie the rows)
#' @return something
#' @author Mark Cowley
#' @exportMethod sort
#' @rdname sort-methods
#' @docType methods
#'
#' @rdname sort-methods
#' @aliases sort,matrix-method
setMethod(
	"sort",
	signature=signature("matrix"),
	function(x, decreasing=FALSE, na.last=NA, FUN, ...) {
		FUN <- match.fun(FUN)
	    
		val <- apply(x, 1, FUN, ...)
		res <- x[order(val, decreasing=decreasing, na.last=na.last), ]

		return(res)
		
	}
)

#' @rdname sort-methods
#' @aliases sort,data.frame-method
setMethod(
	"sort",
	signature=signature("data.frame"),
	function(x, decreasing=FALSE, na.last=NA, FUN, ...) {
		FUN <- match.fun(FUN)
	    
		cols <- colclasses(x) == "numeric"
		res <- x
		res[,cols] <- sort(as.matrix(x[,cols]), decreasing=decreasing, na.last=na.last, FUN=FUN, ...)

		return(res)
		
	}
)

#' @rdname sort-methods
#' @aliases sort,ExpressionSet-method
#' @importClassesFrom Biobase ExpressionSet
setMethod(
	"sort",
	signature=signature("ExpressionSet"),
	function(x, decreasing=FALSE, na.last=NA, FUN, ...) {
		FUN <- match.fun(FUN)
		rn <- rownames(sort(exprs(x), decreasing=decreasing, na.last=na.last, FUN=FUN, ...))
		x <- x[rn, ]
		
		x
	}
)

#' @rdname sort-methods
#' @aliases sort,LumiBatch-method
#' @importClassesFrom lumi LumiBatch
setMethod(
	"sort",
	signature=signature("LumiBatch"),
	function(x, decreasing=FALSE, na.last=NA, FUN, ...) {
		FUN <- match.fun(FUN)
		rn <- rownames(sort(exprs(x), decreasing=decreasing, na.last=na.last, FUN=FUN, ...))
		x <- x[rn, ]
		
		x
	}
)


# #' This is shelved until I have the interface for GCT worked out
# #' and its probable that the ExpressionSet method may just work.
# #' @rdname sort-methods
# #' @aliases sort,GCT-method
# #' @importClassesFrom lumi GCT
# setMethod(
# 	"sort",
# 	signature=signature("GCT"),
# 	function(x, decreasing=FALSE, na.last=NA, FUN, ...) {
# 		FUN <- match.fun(FUN)
# 		rn <- rownames(sort(exprs(x), decreasing=decreasing, na.last=na.last, FUN=FUN, ...))
# 		x[rn, ]
# 	}
# )
# 
# #' Function to reorder the rows of a GCT object
# #' 
# #' take the data in a gct file & reorder the rows by different criteria:\cr
# #' - var: sort by variance, from high to low\cr
# #' - mean: sort by average abundance, from high to low\cr
# #' - median: sort by median abundance, from high to low\cr
# #' - max: sort by max abundance, from high to low\cr
# #' - sum: sort by sum, from high to low (useful if the GCT is an adjacency
# #' matrix/matrix of counts, or scores)\cr
# #' 
# #' @param gct a gct object
# #' @param method see details. Default = 'var'
# #' @param reverse if \code{TRUE}, then reverse the default sort order (see details for
# #'   default sort order)
# #' @author Mark Cowley, 2011-03-16
# #' @export
# #' @seealso \code{\link{reorder.gct.file}}
# reorder.gct <- function(gct, method=c("var", "mean", "median", "max", "sum"), reverse=FALSE) {
# 	N <- ncol(gct)-2
# 	
# 	method <- method[1]
# 	gct <- switch(method,
# 		var={
# 			x <- apply(gct[,3:(N+2)], 1, var)
# 			gct[order(x, decreasing=!reverse), ]
# 		},
# 		mean={
# 			x <- apply(gct[,3:(N+2)], 1, mean)
# 			gct[order(x, decreasing=!reverse), ]
# 		},
# 		median={
# 			x <- apply(gct[,3:(N+2)], 1, median)
# 			gct[order(x, decreasing=!reverse), ]
# 		},
# 		sum={
# 			x <- apply(gct[,3:(N+2)], 1, sum)
# 			gct[order(x, decreasing=!reverse), ]
# 		},
# 		max={
# 			x <- apply(gct[,3:(N+2)], 1, max)
# 			gct[order(x, decreasing=!reverse), ]
# 		},
# 		stop("Unsupported method")
# 	)
# 	
# 	return(gct)
# }
# # CHANGELOG
# # 2012-03-19: reverse parameter did nothing. I've made it actually work.
