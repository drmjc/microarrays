################################################################################
################################################################################
#' replace the controlData slot in a LumiBatch object
#'
#' @param object a LumiBatch object
#' @param value a a \code{data.frame} with first two columns as
#'  \dQuote{controlType} and \dQuote{ProbeID}. The rest columns are the control
#'  probe expression amplitudes for individual samples. To delete the slot,
#'  set \code{value=NULL}.
#' @return something
#' @author Mark Cowley
#' @export
#' @rdname controlDataSetter-methods
#' @docType methods
setGeneric(
	"controlData<-",
	function(object, value) {
		standardGeneric("controlData<-")
	}
)

#' @rdname controlDataSetter-methods
#' @aliases controlData<-,LumiBatch,data.frame-method
#' @importClassesFrom lumi LumiBatch
setMethod(
	"controlData<-",
	signature=signature("LumiBatch", "data.frame"),
	function(object, value) {
		stopifnot(
			ncol(value) == (ncol(object) + 2), 
			identical(colnames(value)[1:2], c("controlType", "ProbeID")),
			identical(colnames(value)[-c(1,2)], sampleNames(object))
		)
		
		object@controlData <- value
		object
	}
)

#' @rdname controlDataSetter-methods
#' @aliases controlData<-,LumiBatch,NULL-method
#' @importClassesFrom lumi LumiBatch
setMethod(
	"controlData<-",
	signature=signature("LumiBatch", "NULL"),
	function(object, value) {
		object@controlData <- NULL
		object
	}
)

################################################################################
################################################################################

# overwrite the method that came from lumi
removeMethod("sampleNames<-", signature=signature("LumiBatch", "ANY"))

#' set sampleNames in LumiBatch objects
#' 
#' set sampleNames in LumiBatch objects, including the QC, vstParameter, transformFun
#'  and controlData slots
#' 
#' @note code came from lumi_2.8.0 & modified by MJC to set the controlData slot
#'  properly
#' 
#' @param object an \code{LumiBatch} object
#' @param value a character vector of column names
#' @return a \code{LumiBatch} object
#' @rdname sampleNamesSetter-methods
#' @aliases sampleNames<-,LumiBatch-method
#' @exportMethod "sampleNames<-"
#' @importMethodsFrom lumi "sampleNames<-"
#' @importMethodsFrom Biobase "sampleNames<-"
#' @importClassesFrom lumi LumiBatch
setMethod("sampleNames<-",
	signature=signature("LumiBatch", "ANY"),
	function(object, value) {
		object <- callNextMethod()
		ddim <- dim(object)
		if (!is.null(object@QC)) {
			QC <- object@QC
			if (!is.null(QC$sampleSummary)) 
				if (ncol(QC$sampleSummary) == ddim[2]) 
					colnames(QC$sampleSummary) <- value
			if (!is.null(QC$BeadStudioSummary)) 
				if (nrow(QC$BeadStudioSummary) == ddim[2]) 
					rownames(QC$BeadStudioSummary) <- value
			object@QC <- QC
		}
		if (!is.null(attr(object, "vstParameter"))) {
			vstParameter <- attr(object, "vstParameter")
			if (!is.null(nrow(vstParameter))) {
				if (nrow(vstParameter) == ddim[2]) {
					rownames(vstParameter) <- value
					transformFun <- attr(object, "transformFun")
					names(transformFun) <- value
					attr(object, "vstParameter") <- vstParameter
					attr(object, "transformFun") <- transformFun
				}
			}
		}
		if (nrow(object@controlData) > 0) {
			# modified by MJC
			colnames(controlData(object)) <- c("controlType", "ProbeID", value)
		}
		return(object)
	}
)

################################################################################
################################################################################

#' subset a LumiBatch object
#' 
#' This is preferable to subset.eSet, since there are extra slots in a 
#' LumiBatch
#' 
#' The support for subsetting the controlData is still inadequately tested.
#' For instance, do all controlData slots have the first 2 columns being
#' "?" and "ProbeID"?
#' 
#' @param x an LumiBatch
#' @param subset logical expression indicating elements or rows to keep:
#'          missing values are taken as \code{FALSE}.
#' @param select logical expression, indicating columns to select from a \code{data.frame}.
#' @param \dots ignored
#' @return an LumiBatch with fewer rows and or columns
#' @author Mark Cowley, 2011-09-01
#' 
#' @importClassesFrom lumi LumiBatch
#' @method subset LumiBatch
#' @S3method subset LumiBatch
#' @rdname subset-LumiBatch-methods
subset.LumiBatch <- function(x, subset, select, ...) {
	if( missing(subset) ) subset <- rep(TRUE, nrow(x))
	if( missing(select) ) select <- rep(TRUE, ncol(x))

	# The [ and ] operators work well for eSet's already.
	res <- x[which(subset), which(select)]
	if(!is.null(res@controlData)) {
		res@controlData <- subset(x@controlData,select=c(TRUE, TRUE, select))
	}
	
	res
}


#' @rdname subset-LumiBatch-methods
#' @aliases subset,LumiBatch-method
#' @export
setMethod("subset",
	signature=signature("LumiBatch"),
	function(x, subset, select, ...) {
		if( missing(subset) ) subset <- rep(TRUE, nrow(x))
		if( missing(select) ) select <- rep(TRUE, ncol(x))

		# The [ and ] operators work well for eSet's already.
		x <- x[which(subset), which(select)]
		
		if (!is.null(x@QC)) {
			x@QC <- subset(x@QC, select=select)
		}
		if (!is.null(attr(x, "vstParameter"))) {
			vstParameter <- attr(x, "vstParameter")
			vstParameter <- subset(vstParameter, subset, select)
			attr(x, "vstParameter") <- vstParameter
		}
		if (!is.null(attr(x, "transformFun"))) {
			transformFun <- attr(x, "transformFun")
			transformFun <- subset(transformFun, subset, select)
			attr(x, "transformFun") <- transformFun
		}
		# NB: different rows to exprs or detection
		if (!is.null(x@controlData)) {
			x@controlData <- subset(x@controlData, select=c(TRUE, TRUE, select))
		}
	}
)

################################################################################
################################################################################

