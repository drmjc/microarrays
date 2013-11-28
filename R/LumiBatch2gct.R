#' convert a LumiBatch to a GCT object
#'
#' convert a LumiBatch to a GenePattern GCT object. There must be an
#' \code{fData(x)$Description} column.
#' 
#' @section TODO:
#' - have a \code{GCT} class, and an \code{as.gct} method to convert the various types to \code{GCT-class}
#' - LumiBatch2res, which would need a LumiBatch2calls method
#' 
#' @param x a LumiBatch object
#' @param description.column the description column within \code{fData(x)} to grab
#'  the data for the resultant Description column in the GCT.
#' 
#' @return a GCT object. this is currently just a \code{data.frame} representation.
#' 
#' @author Mark Cowley, 2012-05-02
#' @export
#' @importClassesFrom lumi LumiBatch
#' @importFrom Biobase fData featureNames exprs
#' 
LumiBatch2gct <- function(x, description.column = "Description") {
    inherits(x, "LumiBatch") || stop("x must be a LumiBatch-class")
    description.column %in% colnames(fData(x)) || stop(sprintf("Could not find description column called '%s' in the featureData slot of x", description.column))
    
    res <- data.frame(
        Name=featureNames(x),
        Description=fData(x)$Description,
        exprs(x)
    )
    res
}

## See DetectionCalls.R
# LumiBatch2res <- function(x, description.column = "Description") {
#     inherits(x, "LumiBatch") || stop("x must be a LumiBatch-class")
#     description.column %in% colnames(fData(x)) || stop(sprintf("Could not find description column called '%s' in the featureData slot of x", description.column))
#     
#     res <- data.frame(
#         Name=featureNames(x),
#         Description=fData(x)$Description
#     )
# 	x <- exprs(x)
# 	p <- DetectionCalls(x)
# 	collate.data.frame(x, p)
#     res
# }
# 
# collate.data.frame <- function(x,y) {
# 	stopifnot(ncol(x)==ncol(y))
# 	stopifnot(nrow(x)==nrow(y))
# 	res <- cbind(x,y)
# 	o <- c(seq(1, ncol(res), 2), seq(2,ncol(res),2))
# 	o <- order(o)
# 	res <- res[,o]
# 	res
# }
# 
# #' Convert DABG P-values into P/M/A calls
# #' 
# #' Convert DABG P-values into P/M/A calls, by hard-thresholding at two
# #' threshold, one for "Present", and "Marginal"
# #' 
# #' @param dabg a data.frame of dabg calls. see ?import.APT
# #' @param Pthresh the p-values at which to threshold the DABG p-values
# #' @param Mthresh the p-values at which to threshold the DABG p-values
# #' @return a data.frame, with same dim, and dimname as dabg, but with values of
# #'   "P", "M", "A"
# #' @author Mark Cowley, 2008-05-26
# #' @export
# dabg2calls <- function(dabg, Pthresh=1e-05, Mthresh=0.001) {
# 	calls <- as.data.frame(matrix("A", nrow(dabg), ncol(dabg)), stringsAsFactors=FALSE)
# 	dimnames(calls) <- dimnames(dabg)
# 	calls[dabg < Pthresh] <- "P"
# 	calls[dabg < Mthresh & calls != "P"] <- "M"
# 	
# 	calls
# }
# 
