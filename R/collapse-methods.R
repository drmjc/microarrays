#' collapse
#' 
#' collapse a dataset to have just 1 row per unique value of key.
#'
#' @details given an object with data, and 1 'key' column, containing possibly non-unique
#' identifiers, create a result, which has unique values for those keys. The 'key' column
#' could be a gene symbol, and each row could be a probe; the goal is to convert a 
#' 1-row-per-probe table into 1-row-per-gene.
#' 
#' @param x an object
#' @param decreasing logical: should the sort order be increasing or decreasing?
#' @param na.last for controlling the treatment of \code{NA}s.  If \code{TRUE}, missing
#'   values in the data are put last; if \code{FALSE}, they are put
#'   first; if \code{NA}, they are removed.
#' @param FUN a function used to determine which row to include, when there are multiple 
#'   rows with the same 'column' value. eg \code{mean, min, max, var, sd, mad, prod}
#' @param column the column name within \code{x} which contains the keys. see details
#' @param \dots arguments passed to FUN
#' 
#' @return something
#' @author Mark Cowley
#' @exportMethod collapse
#' @rdname collapse-methods
#' @docType methods
setGeneric(
	"collapse",
	function(x, decreasing=TRUE, na.last=NA, FUN=mean, column, ...) {
		standardGeneric("collapse")
	}
)

#' @rdname collapse-methods
#' @aliases collapse,LumiBatch,missing,missing,missing,missing-method
setMethod(
	"collapse",
	signature=signature("LumiBatch", "missing", "missing", "missing", "missing"),
	function(x, decreasing=TRUE, na.last=NA, FUN=mean, column="SymbolReannotated", ...) {
		collapse(x, decreasing=TRUE, na.last=NA, FUN=mean, column="SymbolReannotated")
	}
)

#' @section collapse,LumiBatch and column:
#' The \code{column} should be the name of a column found in the \code{fData(x)} slot.
#' hint: \code{colnames(fData(x)), or fvarLabels(x)}
#' 
#' @section collapse,LumiBatch and Missing values:
#' Often some probes don't have a genesymbol, and should thus have an \code{NA} in this
#' column. \code{na.list} controls what to do with the gene-less probes. 
#' \code{na.last=TRUE} keeps the probes & moves them to the bottom of the result. 
#' \code{na.last=FALSE} keeps the probes & leaves them in the sorted order specified by \code{FUN}.
#' Note this is different to other implementations of \code{na.last}.
#' \code{na.last=NA (default)} discards these probes from the result.
#' 
#' @rdname collapse-methods
#' @aliases collapse,LumiBatch,logical,logical,function,character-method
#' @importFrom plyr join
#' 
#' @examples
#' \dontrun{
#'  load("Rmisc/x.averaged.RDa.gz")
#' 	tmp <- collapse(x.averaged, FUN=var, decreasing=TRUE, na.last=FALSE, "SymbolReannotated")
#' }
#' 
#' require(lumi)
#' data(example.lumi)
#' featureData(example.lumi)$GeneSymbol <- c("TP53", "BRCA1", "BRCA2", "KRAS", "SMAD4")
#' collapse(example.lumi, TRUE, NA, mean, "GeneSymbol")
#' 
setMethod(
	"collapse",
	signature=signature("LumiBatch", "logical", "logical", "function", "character"),
	function(x, decreasing=TRUE, na.last=NA, FUN=mean, column="SymbolReannotated", ...) {
		FUN <- match.fun(FUN)
		if( ! column %in% fvarLabels(x) ) {
			stop(sprintf("Can't find Gene Symbol column called %s, within fData(x). available columns are:\n%s", column, paste(fvarLabels(x), collapse=", ")))
		}
		fData(x)$featureNames <- featureNames(x)
		
		x <- sort(x, decreasing=decreasing, na.last=na.last, FUN=FUN, ...)
		genes <- na.omit(unique(fData(x)[,column]))
		genes.idx <- match(genes, fData(x)[,column])
		
		#
		# remember the other possibilities for each probe-to-gene mapping
		# 
		tmp <- data.frame(fData(x)[,column], fData(x)$featureNames, stringsAsFactors=FALSE)
		colnames(tmp) <- c(column, "PossibleProbes")
		b <- collapse.rows(tmp, 1, 2, ", ")
		fData(x) <- join(fData(x), b, column, "left")
		fData(x)$PossibleProbes[is.na(fData(x)$PossibleProbes)] <- fData(x)$featureNames[is.na(fData(x)$PossibleProbes)]
		fData(x)$featureNames <- NULL
		
		#
		# handle the NA's, ie probes with no gene symbol.
		# 
		if( is.na(na.last) ) { # NOP
			na.last <- NA
			new.featureNames <- fData(x)[,column]
		}
		else if( na.last ) {
			# move NA probes down to the bottom of the table
			na.idx <- which(is.na(fData(x)[,column]))
			new.featureNames <- c(fData(x)[genes.idx, column], featureNames(x)[na.idx])
			genes.idx <- c(genes.idx, na.idx)
		}
		else if( !na.last ) {
			# leave NA probes in-place & discard the duplicates.
			genes.idx <- which(!duplicated(fData(x)[,column]) | is.na(fData(x)[,column]))
			new.featureNames <- ifelse(is.na(fData(x)[genes.idx,column]), featureNames(x)[genes.idx], fData(x)[genes.idx,column])
		}
		
		res <- x[genes.idx, ]
		featureNames(res) <- new.featureNames

		return( res )
	}
)

# #' Collapse a GCT object
# #' 
# #' It's common for microarrays to have multiple probes per gene. They tend to
# #' represent different isoforms.
# #' Most geneset testing is done at the gene symbol level & ignores isoforms,
# #' so you need to choose 1 probe
# #' for each gene.
# #' How?
# #' 2 common approaches are to take the most abundant probe, or the most
# #' variable probe, considered
# #' across the cohort.\cr
# #' I quite like doing t-stats on each gene & selecting the best performing
# #' probe - ie the one with the largest t-stat in
# #' either direction. Why? On the Affy 133+2 array, there can be lots of poor
# #' probes for each gene. If 5 probes for a gene
# #' have these t-stats: 1.2, 0.9, 0.1, -0.1, -10; then IMO, the one that scored
# #' -10 is the best probe, since it had a really
# #' strong t-stat score. thus method="maxabs" combined with a rnk.file (NB NOT
# #' implemented right now.)
# #' 
# #' @param gct a GCT object
# #' @param chip a CHIP object
# #' @param rnk [optional] path to a rnk file (eg a t-statistic for each probe,
# #'   where you want to select best probe from this score) NB currently UNUSED
# #' @param method \dQuote{mean}, \dQuote{median} select the probe with highest 
# #'   average/median level, or \dQuote{var}: select the probe with highest variance 
# #'   across samples; \dQuote{maxabs} select the probe with the large absolute score in the rnk
# #'   file (see details).
# #' @param reverse [default=FALSE] reverse the ordering selected by method arg.
# #'   so instead of most variable, it would be least variable.
# #' @param filter Filter out (ie exclude) those probes that don't have a gene
# #'   symbol (as determined by probes that have a gene symbol of \code{NA}, \dQuote{},
# #'   \dQuote{---}, or \dQuote{NA}.)
# #' 
# #' @return A gct object with 1 row per gene symbol & now the \sQuote{probe ids} in
# #'   column 1 are actually gene symbols.
# #' @author Mark Cowley, 2011-02-27
# #' @export
# #' @rdname collapse-methods
# #' @aliases collapse,GCT-method
# setMethod(
# 	"collapse",
# 	signature=signature("GCT"),
# 	function(x, method=c("var", "mean", "median", "max"), reverse=FALSE, chip, rnk=NULL, filter=FALSE, ...) {
# 		gct <- x
# 		method <- match.arg(method)
# 		!missing(gct) && all(c("Name", "Description") %in% colnames(gct)) || stop("gct should be a gct object")
# 		!missing(chip) || stop("chip must be specified")
# 
# 		if( !all(gct[,1] %in% chip[,1]) ) {
# 			stop("Some ID's in gct are not in the chip file.")
# 		}
# 		else {
# 			chip <- chip[match(gct[,1], chip[,1]), ]
# 		}
# 
# 		#
# 		# reorder the rows of the gct so that better probes are higher up that worse probes.
# 		# (better is set according to the given 'method')
# 		#
# 		gct <- reorder.gct(gct, method=method, reverse=reverse)
# 		chip <- chip[match(gct[,1], chip[,1]), ]
# 		gct$ORDER <- 1:nrow(gct) # used to integrate the rows that do and do NOT have probe symbols.
# 
# 		# split the data into those probes that do & do not have Gene Symbols
# 		probes.no.sym   <- chip[,2] == "NA" | chip[,2] == "---" | chip[,2] == "" | is.na(chip[,2])
# 		gct.hasSymbols  <- gct[!probes.no.sym, ]
# 		gct.noSymbols   <- gct[probes.no.sym, ]
# 		chip.hasSymbols <- chip[!probes.no.sym, ]
# 		chip.noSymbols  <- chip[probes.no.sym, ]
# 
# 		# for the probes that DO have a valid gene symbol, choose the best probe.
# 		uSymbols <- unique(chip.hasSymbols[,2])
# 		m <- match(uSymbols, chip.hasSymbols[,2])
# 		res <- data.frame( chip.hasSymbols[m, 2:3], 
# 						   gct.hasSymbols[m, 3:ncol(gct.hasSymbols)], 
# 						   check.names=FALSE )
# 		colnames(res)[1:2] <- colnames(gct)[1:2]
# 		rownames(res) <- res[,1]
# 		res[,2] <- sprintf("%s (Probe:%s)", res[,2], gct.hasSymbols[m,1])
# 
# 		# for the probes that DO NOT have a valid gene symbol, use the same information.
# 		if( !filter && (nrow(gct.noSymbols) > 0) ) {
# 			# Reannotate the first 2 columns using the chip file
# 			gct.noSymbols[,2] <- chip2description(chip=chip.noSymbols)
# 			# gct.noSymbols[,2] <- sprintf("%s (Probe:%s)", gct.noSymbols[,2], gct.noSymbols[,1])
# 			res <- rbind(res, gct.noSymbols)
# 			res <- res[order(res$ORDER), ]
# 		}
# 
# 		res$ORDER <- NULL
# 
# 		res
# 	}
# )
# # CHANGELOG
# # 2012-03-19: added (Probe:xyz) into the Description column.
# # 2012-07-05: migrated from collapse.gct
