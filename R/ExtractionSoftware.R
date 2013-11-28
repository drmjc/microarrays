#' ExtractionSoftware
#'
#' Determine which extraction software was used to create the 
#' Illumina TXT file. The first header of the TXT file looks something like this:
#' \sQuote{Illumina Inc. GenomeStudio version 1.8.0}, or 
#' \sQuote{Illumina Inc. BeadStudio version 1.4.0.1}. This information is embedded 
#' within the LumiBatch object.
#' 
#' @param x a LumiBatch object
#' 
#' @return \sQuote{GenomeStudio}, or \sQuote{BeadStudio}
#' 
#' @author Mark Cowley
#' @exportMethod ExtractionSoftware
#' @rdname ExtractionSoftware-methods
#' @docType methods
#' @examples
#' require(lumi)
#' data(example.lumi)
#' ExtractionSoftware(example.lumi)
setGeneric(
	"ExtractionSoftware",
	function(x) {
		standardGeneric("ExtractionSoftware")
	}
)

#' @rdname ExtractionSoftware-methods
#' @aliases ExtractionSoftware,LumiBatch-method
setMethod(
	"ExtractionSoftware",
	signature=signature("LumiBatch"),
	function(x) {
		hdr <- x@experimentData@other[[1]][1]
		if( length(grep("GenomeStudio", hdr))>0 )
			res <- "GenomeStudio"
		else if( length(grep("BeadStudio", hdr)>0) )
			res <- "BeadStudio"
		else {
			msg <- c("Couldn't determine extraction software. Here's the header:", x@experimentData@other[[1]])
			msg <- paste(msg, collapse="\n")
			stop(msg)
		}
		res
	}
)
