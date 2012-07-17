#' add a lowess curve to a plot created by plotMA 
#' 
#' @param MA a \code{list}, \code{MAList}, \code{RGList}, \code{MArrayLM}
#' @param array an index into the MA object
#' @param col colur
#' @param \dots arguments passed to \code{\link[stats]{lowess}}
#' @author Mark Cowley, 2008-08-19
#' @export
#' @importFrom limma MA.RG
#' @importClassesFrom limma MAList RGList MArrayLM
#' @importFrom stats lowess
#' 
#' @examples
#' \dontrun{
#' plotMA(RG, 1)
#' lowessMA(RG, 1, col="red")
#' }
lowessMA <- function(MA, array=1, col="purple", ...) {
    # stolen from plotMA
    if (class(MA) == "list") 
        MA <- new("MAList", MA)
    if (is(MA, "RGList")) {
        MA <- MA.RG(MA[, array])
        array <- 1
    }
    if (is(MA, "MAList")) {
        x <- as.matrix(MA$A)[, array]
        y <- as.matrix(MA$M)[, array]
    }
    else if (is(MA, "MArrayLM")) {
        if (is.null(MA$Amean)) 
            stop("MA-plot not possible because Amean component is absent.")
        x <- MA$Amean
        y <- as.matrix(MA$coef)[, array]
    }
    else {
        MA <- as.matrix(MA)
        narrays <- ncol(MA)
        if (narrays < 2) 
            stop("Need at least two arrays")
        if (narrays > 5) 
            x <- apply(MA, 1, median, na.rm = TRUE)
        else x <- rowMeans(MA, na.rm = TRUE)
        y <- MA[, array] - x
    }

	nas <- union( which(is.na(x)), which(is.na(y)) )
	ok <- setdiff(1:length(x), nas)
    lines( lowess(x[ok],y[ok]), col=col, ... )
}
