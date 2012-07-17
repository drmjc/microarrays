#' From a toptable which contains probe ID's in the ID column, determine the
#' best probe for each gene
#' 
#' If multiple probes map to a gene symbol, then just 1 row will be returned -
#' the one with the max(abs(tstat))
#' If the probe has no symbol (or is "---", "", NA, NULL, symbols.ignore.list)
#' then it will also be exported.
#' 
#' @param tt a toptable
#' @param probe2gene a 2 column data frame with probes and gene symbols
#' @param toupper if TRUE, convert the gene symbol to UPPERCASE; FALSE leaves
#'   it untouched
#' @param symbols.ignore.list Probes with these values as gene symbols, as well
#'   as NA and NULL will not be merged to the same value, otherwise the 1000
#'   probes that map to "---" will be collapsed to a single 'gene'
#' @param only.genes After collapsing multiple probes to one, if TRUE, then
#'   only the probes with real gene symbols will be exported; if FALSE, all
#'   collapsed probes with or without a gene symbol will be exported.
#' @author Mark Cowley, 2008-12-08
#' @export
calc.best.probe.topTable <- function(tt, probe2gene, toupper=TRUE, symbols.ignore.list=c("---", ""), only.genes=FALSE) {
	# clean up the probe 2 gene symbols section.
	idx <- is.na(probe2gene[,2]) | is.null(probe2gene[,2]) | nchar(probe2gene[,2])==0
	if( any(idx) ) {
		probe2gene[idx,2] <- "---"
	}
	if( all(probe2gene[,2] %in% symbols.ignore.list) ) {
		res <- data.frame(ProbeSetID=tt$ID, Gene.Symbol="", Nprobes=1, Possible.Probes=tt$ID, stringsAsFactors=FALSE)
		return(res)
	}

	colnames(probe2gene) <- c("ID", "Gene.Symbol")
	if( any(colclasses(probe2gene) != "character") ) colclasses(probe2gene) <- rep("character", ncol(probe2gene))
	if( toupper )
		probe2gene[,2] <- toupper(probe2gene[,2])
	
	# re-order the probe 2 symbol mapping to match the tt order.
	probe2gene <- probe2gene[match(tt$ID, probe2gene[,1]), ]
	
	# trim out those probes with no gene symbol
	noSymbols    <- probe2gene[  probe2gene[,2] %in% symbols.ignore.list, ]
	probe2gene <- probe2gene[! probe2gene[,2] %in% symbols.ignore.list, ]

	# collapse those probes that have a gene symbol
	sym2indices <- as.list(vector2hashtable(probe2gene[,2], warn=FALSE))
	first <- sapply(sym2indices, "[", 1)
	res <- probe2gene[first, ]
	res$Nprobes <- sapply(sym2indices, length)
	res$Possible.Probes <- sapply(sym2indices, function(x) {
		paste(probe2gene[x,1], collapse=", ")
	})
	
	# add back the probes without gene symbols.
	if( !only.genes && nrow(noSymbols) > 0 ) {
		noSymbols$Nprobes <- 1
		noSymbols$Possible.Probes <- noSymbols$ID
		res <- rbind(res, noSymbols)
		# reorder the rows of res to match the original order of rows from tt, allowing for the loss of probes due to collapsing.
		res <- res[match(res$ID, tt$ID[tt$ID %in% res$ID]), ]
	}
	res
}
