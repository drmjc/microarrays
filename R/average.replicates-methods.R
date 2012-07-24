#' Average replicate columns
#' 
#' @description Average the replicate samples within a dataset \code{x}, based on the groupings indicated
#' in by \code{classes}, where identical values indicate samples to be \code{averaged}.
#' Note that \code{averaged} is a bit of a generalisation, since some dataslots, eg
#' \code{detection} within \code{LumiBatch} objects contains p-values, which are not averaged per se,
#' see the \code{average.pval.replicates} S4 methods.
#' 
#' @details \code{x} can be either a \code{matrix}, \code{data.frame}, \code{eSet}, \code{LumiBatch},
#' \code{AnnotatedDataFrame}. Note that for most of these data types, each sample is assumed to be
#' in a column. \code{AnnotatedDataFrame} has 1 row per sample, thus the result will have fewer rows
#' than \code{x}.
#' 
#' \code{classes} can be a character or factor vector. Note that the \code{sampleNames}
#' of the result are derived from the unique values, or levels in the \code{classes} vector
#'
#' @param x a \code{matrix}, \code{data.frame}, \code{eSet}, \code{LumiBatch} or \code{AnnotatedDataFrame} of data
#' 
#' @param classes a \code{character} or \code{factor}: the sample classes, 1 per sample in \code{x}, with as many
#'   unique values or \code{levels} as there are unique samples.
#' 
#' @return The same data type as \code{x}, but with fewer samples, where the replicates have been averaged.
#' 
#' @author Mark Cowley, 2011-09-01
#' 
#' @exportMethod average.replicates
#' @rdname average.replicates-methods
#' @docType methods
#' 
#' @examples
#' \dontrun{
#' load("Rmisc/x.norm.RDa.gz")
#' classes <- sub("\\.[12]$", "", sampleNames(x.norm))
#' 
#' # Matrix example
#' mat <- exprs(x.norm)
#' avg <- average.replicates(mat, classes)
#' 
#' # data.frame example
#' df <- as.data.frame(mat)
#' df.avg <- average.replicates(df, classes)
#' 
#' # LumiBatch example
#' x.av <- average.replicates(x.norm, classes)
#' 
#' # AnnotatedDataFrame example
#' adf <- PhenoData(x.norm)
#' adf.avg <- average.replicates(adf, classes)
#' adf.avg
#' as(adf.avg, "data.frame")
#' }
setGeneric(
	"average.replicates",
	function(x, classes) {
		standardGeneric("average.replicates")
	}
)

#' @details \code{matrix} will return a \code{matrix}, with 1 column per 
#' \code{levels(classes)}, in the
#' same order as \code{classes}. Note this is different
#' 
#' @rdname average.replicates-methods
#' @aliases average.replicates,matrix,factor-method
setMethod(
	"average.replicates",
	signature=signature("matrix", "factor"),
	function(x, classes) {
		length(classes) == ncol(x) || stop("length(classes) != ncol(x)")
		design <- model.matrix(~0+classes)
		colnames(design) <- levels(classes)
		
		fit <- lm.fit(design, t(x))
		res <- t(fit$coefficients)
		res
	}
)

#' @rdname average.replicates-methods
#' @aliases average.replicates,matrix,character-method
setMethod(
	"average.replicates",
	signature=signature("matrix", "character"),
	function(x, classes) {
		classes <- factor(classes, levels=unique(classes))
		average.replicates(x, classes)
	}
)

#' @rdname average.replicates-methods
#' @aliases average.replicates,data.frame,ANY-method
setMethod(
	"average.replicates",
	signature=signature("data.frame", "ANY"),
	function(x, classes) {
		res <- average.replicates(as.matrix(x), classes)
		res <- as.data.frame(res, stringsAsFactors=FALSE)
		res
	}
)

#' @rdname average.replicates-methods
#' @aliases average.replicates,ANY,missing-method
setMethod(
	"average.replicates",
	signature=signature("ANY", "missing"),
	function(x, classes) {
		stop("You must supply a classes vector.")
	}
)
################################################################################

#' @details \code{eSet} is a virtual object, thus all sub-classes of eSet, like
#' ExpressionSet or AffyBatch should end up calling this method. This method
#' will work as long as there are no additional slots defined in the class
#' (which is why I have a separate LumiBatch-specific method). \code{average.replicates}
#' on \code{eSet}'s does: (1) resizes \emph{all} elements within \code{AssayData},
#' (which should all be matrix-like, and includes exprs, se.exprs & any other items found), 
#' as well as (2) \code{phenoData} and (3) \code{protocolData}, which indeed are the only 
#' elements defined in the \code{eSet} interface, which depend on the number of samples.\cr
#' IF you get an error indicating that not all slots in the result contain the correct number
#' of samples, then you may need to create another S4 method for that particular \code{eSet-sub-class}. Use
#' the method implementation for \code{eSet} as a starting point, and then the method for \code{LumiBatch} 
#' as a worked example of an \code{eSet-sub-class} with a few additional slots.
#' For (2,3) above, which are AnnotatedDataFrame's the replicate samples are averaged using the 
#' \code{average.replicates} for AnnotatedDataFrame method.
#' 
#' @rdname average.replicates-methods
#' @importClassesFrom Biobase eSet
#' @importMethodsFrom Biobase assayData phenoData "phenoData<-" protocolData "protocolData<-"
#' @importFrom Biobase assayDataElementNames assayDataElementReplace
#' @aliases average.replicates,eSet,factor-method
setMethod(
	"average.replicates",
	signature=signature("eSet", "factor"),
	function(x, classes) {
		length(classes) == ncol(x) || stop("length(classes) != ncol(x)")
		
		res <- x
		assay.dat.elements <- assayDataElementNames(x)
		for(i in seq(along=assay.dat.elements)) {
			assay.dat.elem.name <- assay.dat.elements[i]
			assay.dat.elem <- average.replicates(assayData(res)[[assay.dat.elem.name]], classes)
			res <- assayDataElementReplace(res, assay.dat.elem.name, assay.dat.elem)
		}
		
		pheno.dat <- average.replicates(phenoData(x), classes)
		phenoData(res) <- pheno.dat

		proto.dat <- average.replicates(protocolData(x), classes)
		protocolData(res) <- proto.dat
		
		# force a re-evaluation of the object
		sampleNames(res) <- levels(classes)
		
		res
	}
)

#' @rdname average.replicates-methods
#' @importClassesFrom Biobase eSet
#' @aliases average.replicates,eSet,character-method
setMethod(
	"average.replicates",
	signature=signature("eSet", "character"),
	function(x, classes) {
		classes <- factor(classes, levels=unique(classes))
		average.replicates(x, classes)
	}
)

################################################################################
###################### LumiBatch ###############################################
################################################################################

#' @details \code{\linkS4class{LumiBatch}} objects contain a few additional slots in addition
#' to those required by the \code{\linkS4class{eSet}} interface: \code{QC}, \code{controlData}.
#' Also, one of the 
#' \code{assayData} slots, \code{detection} contains detection p-values which
#' shouldn't really just be averaged.
#' This uses \code{average.pval.replicates} to properly average the detection pvalues
#' from \code{\linkS4class{LumiBatch}} objects with replicates.
#' 
#' @rdname average.replicates-methods
#' 
#' @importClassesFrom Biobase ExpressionSet
#' @importClassesFrom lumi LumiBatch
#' @importMethodsFrom lumi detection "detection<-" getHistory
#' @aliases average.replicates,LumiBatch,factor-method
#' 
#' @seealso \code{\link{average.pval.replicates}}
setMethod(
	"average.replicates",
	signature=signature("LumiBatch", "factor"),
	function(x, classes) {
		sub.time   <- as.character(Sys.time())
		length(classes) == ncol(x) || stop("length(classes) != ncol(x)")
		
		# update all of the ExpressionSet elements
		res <- as(x,"ExpressionSet")
		res <- average.replicates(res, classes)
		res <- as(res, "LumiBatch")
		# res <- x
		
		# and now, all the LumiBatch elements:
		
		# update the QC slot (res@QC$history is updated later)
		if( !is.null(x@QC) ) {
			res@QC <- x@QC
			QC.ss <- average.replicates(res@QC$sampleSummary, classes)
			res@QC$sampleSummary <- QC.ss
		}
		
		# update the controlData slot
		if( !is.null(x@controlData) && ncol(x@controlData) == ncol(x)+2 ) {
			cd <- data.frame(
				x@controlData[,c(1:2)], 
				average.replicates(x@controlData[,-c(1:2)], classes),
				stringsAsFactors=FALSE, check.names=FALSE, check.rows=FALSE,
				row.names=rownames(x@controlData)
			)
			res@controlData <- cd
		}
		
		# were detection pvalues defined?
		if( !is.null(detection(x)) ) {
			detection(res) <- average.pval.replicates(detection(x), classes)
		}
		
		################################################################################
		# update the history & QC history tables
		.update.history <- function(old.history, start.time, end.time, cmd) {
			if (is.null(old.history$lumiVersion)) 
				old.history$lumiVersion <- rep(NA, nrow(old.history))
			lumiVersion <- packageDescription("lumi")$Version
			res <- rbind(
				old.history, 
				data.frame(submitted = start.time, 
						   finished = end.time, 
						   command = cmd, 
						   lumiVersion = lumiVersion
						  )
			)
			res
		}
		
		cmd <- capture.output(print(match.call(average.replicates)))
		finished.time <- as.character(Sys.time())
		res@history   <- .update.history(getHistory(x), sub.time, finished.time, cmd)
		if( !is.null(res@QC) ) {
			res@QC$history  <- .update.history(x@QC$history,  sub.time, finished.time, cmd)
		}
		################################################################################
		
		# this should force a re-evaluation of the res object for correct dimensions.
		sampleNames(res) <- levels(classes)
		
		res
	}
)

#' @rdname average.replicates-methods
#' @importClassesFrom lumi LumiBatch
#' @aliases average.replicates,LumiBatch,character-method
setMethod(
	"average.replicates",
	signature=signature("LumiBatch", "character"),
	function(x, classes) {
		classes <- factor(classes, levels=unique(classes))
		average.replicates(x, classes)
	}
)

################################################################################


################################################################################
#' @details \code{AnnotatedDataFrame} contains 1 row per sample, there's nothing in
#' the spec as to required columns. 
#' The default for \code{LumiBatch} objects is a single
#' column called sampleID, but this does not have to be the case.
#' It's not clear how to average \code{character} values, so we retain
#' the rows for the first sample matching each unique \code{class}.
#' @rdname average.replicates-methods
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @aliases average.replicates,AnnotatedDataFrame,factor-method
setMethod(
	"average.replicates",
	signature=signature("AnnotatedDataFrame", "factor"),
	function(x, classes) {
		data <- as(x, "data.frame")
		final.classes <- levels(classes)
		if(ncol(data) == 0) {
			data <- data.frame(matrix(NA,length(final.classes), 0))
		}
		else {
			data <- data[match(final.classes, as.character(classes)),,drop=FALSE ]
			rownames(data) <- final.classes
		}
		
		res <- x
		pData(res) <- data
		sampleNames(res) <- final.classes
		
		res
	}
)

#' @rdname average.replicates-methods
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @aliases average.replicates,AnnotatedDataFrame,character-method
setMethod(
	"average.replicates",
	signature=signature("AnnotatedDataFrame", "character"),
	function(x, classes) {
		classes <- factor(classes, levels=unique(classes))
		average.replicates(x, classes)
	}
)
################################################################################
