#' DetectionCalls
#' 
#' From microarray data, extract the detection p-values, and convert to P/M/A calls, ie
#' Present/Marginal/Absent.
#' 
#' @section Affymetrix Genome/IVT arrays:
#' Older generation Affymerix arrays use \code{affy::\link{mas5calls}}, where the thresholds are
#' \code{c(0.04, 0.06)}.
#' 
#' @section Affymetrix ST arrays:
#' ST arrays, that use the DABG procedure for detection pvalues need to have quite stringent
#' p-value thresholds. We recommend \code{c(1e-05, 0.001)}.
#'
#' @section LumiBatch objects:
#' As far as I can tell, BeadStudio and GenomeStudio differ in how they report detection p-values.
#' BeadStudio (older) has low values indicating expression; GenomeStudio reports 1-p, thus
#' high values indicate high expression. This is true at least for BeadStudio 1.4.0.1, and
#' GenomeStudio 1.8.0 files. In either case, the values for thresh should be small (default=c(0.01, 0.05)),
#' and for GenomeStudio, the logic will be reversed such that probes with p-val > 0.95 are
#' Marginal, and p-val > 0.99 are Present.
#' 
#' @param x an ExpressionSet, or LumiBatch object
#' @param thresh a numeric(2) where the values represent the p-value
#' thresholds for Present and Marginal, respectively.
#' Some guidance re Illumina detection pvals here \url{http://www.illumina.com/Documents/products/technotes/technote_gene_expression_data_quality_control.pdf}
#' 
#' @return a character matrix, same dimension as x
#' 
#' @author Mark Cowley
#' @exportMethod DetectionCalls
#' @rdname DetectionCalls-methods
#' @docType methods
#' @examples
#' require(lumi)
#' data(example.lumi)
#' head(DetectionCalls(example.lumi))
setGeneric(
	"DetectionCalls",
	function(x, thresh) {
		standardGeneric("DetectionCalls")
	}
)

#' @rdname DetectionCalls-methods
#' @aliases DetectionCalls,ExpressionSet,missing-method
setMethod(
	"DetectionCalls",
	signature=signature("ExpressionSet", "missing"),
	function(x, thresh) {
		DetectionCalls(x, c(1e-05,0.001))
		# thresh <- c(1e-05,0.001)
		# callNextMethod()
	}
)

#' @rdname DetectionCalls-methods
#' @aliases DetectionCalls,ExpressionSet,numeric-method
setMethod(
	"DetectionCalls",
	signature=signature("ExpressionSet", "numeric"),
	function(x, thresh) {
		calls <- matrix("A", nrow(x), ncol(x))
		calls[detection(x) < thresh[2]] <- "M"
		calls[detection(x) < thresh[1]] <- "P"
		dimnames(calls) <- list(featureNames(x), sampleNames(x))
		calls
	}
)


#' @rdname DetectionCalls-methods
#' @aliases DetectionCalls,LumiBatch,missing-method
setMethod(
	"DetectionCalls",
	signature=signature("LumiBatch", "missing"),
	function(x, thresh) {
		DetectionCalls(x, c(0.01, 0.05))
		# thresh=c(0.01, 0.05)
		# callNextMethod()
	}
)

#' @rdname DetectionCalls-methods
#' @aliases DetectionCalls,LumiBatch,numeric-method
setMethod(
	"DetectionCalls",
	signature=signature("LumiBatch", "numeric"),
	function(x, thresh) {
		calls <- matrix("A", nrow(x), ncol(x))
		if( ExtractionSoftware(x) == "GenomeStudio" ) {
			calls[detection(x) > 1-thresh[2]] <- "M"
			calls[detection(x) > 1-thresh[1]] <- "P"
		}
		else if( ExtractionSoftware(x) == "BeadStudio" ) {
			calls[detection(x) < thresh[2]] <- "M"
			calls[detection(x) < thresh[1]] <- "P"
		}
		else {
			stop("Can't determine the ExtractionSoftware used to create the LumiBatch object, so can't determine the detection calls")
		}
		dimnames(calls) <- list(featureNames(x), sampleNames(x))
		calls
	}
)
