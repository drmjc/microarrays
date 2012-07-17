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
#' @return a GCT object. this is currently just a \code{data.frame} representation.
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
