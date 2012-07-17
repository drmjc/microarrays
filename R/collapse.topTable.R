#' Take a topTable with one row per probeset, and collapse it to one row per
#' gene symbol.
#' 
#' Probes either map to Gene Symbols, or not & some Genes can have multiple
#' probes mapping to them (due to array design, or isoforms).
#' This method will choose the best probe per gene, by selecting that with the
#' largest abs(t.stat), since that probe 'performed' the best, (can be either
#' up or down-regulated). In my experience, for those genes with lots of
#' probes, many do not hyb well (eg on the HGU133+2 array, due to poor probe
#' design), and one will hyb better than the others & thus produce larger +/-
#' t.stat. This is the best peforming probe & should be selected to represent
#' that gene. Other alternatives include median probe, max.avg, ... (not yet
#' implemented.)
#' Probes that don't map to symbols can either be filtered out
#' (only.genes=TRUE), or left in (only.genes=FALSE). Probes that
#' don't map to gene symbols are identified via the probe2gene table & either
#' have NA, NULL, or symbols.ignore.list in the 2nd column of the probe2gene.
#' You can add custom values to the symbols.ignore.list, eg c("N/A", "missing",
#' "blank", "neg.con", ...)
#' 
#' @param tt a toptable, or a list of toptables
#' @param probe2gene a 2 column data.frame, 1 row per probe, mapping probe ID's
#'   (column 1) to gene symbols (column 2)
#' @param toupper convert the gene symbol in column 2 to UPPER case?
#' @param symbols.ignore.list a vector of symbols that are used to identify
#'   probes with no gene. NB: the special values NA and NULL are always used,
#'   so no need to specify these.
#' @param only.genes if TRUE, then filter out the probes without a gene symbol.
#'   if FALSE, then probes with or without gene symbols are exported.
#' @param verbose Helpful messages?
#' @return a data.frame like tt, but with some new columns: "Gene.Symbol" in
#'   pos 2, "Nprobes" = the number of possible probes for the given gene,
#'   "Possible.Probes" = a comma separated character(1) of the possible probes
#'   for the given gene. nrows = number of unique genes + if(only.genes=FALSE)
#'   the number of probes with no gene. row ordering is given by the original
#'   tt's row ordering.
#' @author Mark Cowley, 2011-02-22
#' @export
collapse.topTable <- function(tt, probe2gene, toupper=FALSE, symbols.ignore.list=c("---", ""), only.genes=FALSE, verbose=TRUE) {
	if( is.list(tt) && !is.data.frame(tt) ) {
		return(lapply(tt, collapse.topTable, probe2gene=probe2gene, toupper=toupper, symbols.ignore.list=symbols.ignore.list, only.genes=only.genes, verbose=verbose))
	}
	
	n.probes <- nrow(tt)

	if( any(nchar(probe2gene[,2]) > 256) ) {
		idx <- which(nchar(probe2gene[,2]) > 256)
		probe2gene[idx,2] <- paste(str.left(probe2gene[idx,2], 252), "...")
	}
	
	# this will be used to keep the table sorted properly.
	tt$old_order <- 1:nrow(tt)
	
	probe2gene <- probe2gene[probe2gene[,1] %in% tt$ID, ]
	colnames(probe2gene) <- c("ID", "Gene.Symbol")
	
	# deal with probes that map to gene symbols first.
	best.probe <- calc.best.probe.topTable(tt, probe2gene, toupper=toupper, symbols.ignore.list=symbols.ignore.list, only.genes=only.genes)
	tt.genes <- merge(tt, best.probe, by.x="ID", by.y="ID", all.x=FALSE, all.y=TRUE, sort=FALSE)
	tt.genes <- move.column(tt.genes, "Gene.Symbol", 2)

	n.genes <- sum(tt.genes$Gene.Symbol != "---")
	n.extra <- sum(tt.genes$Gene.Symbol == "---")

	res <- tt.genes
	res <- res[order(as.numeric(res$old_order)), ]
	rownames(res) <- 1:nrow(res)
	res$old_order <- NULL

	# n.genes <- nrow(tt.genes)
	# n.extra <- 0
	# if( !only.genes ) {
	# 	idx <- which(probe2gene[,2] %in% symbols.ignore.list)
	# 	if( length(idx) > 0 ) {
	# 		ids <- probe2gene[idx, 1]
	# 		tt.probes <- tt[match(ids, tt$ID), ]
	# 		res <- rbind.smart(tt.genes, tt.probes)
	# 		idx <- which(is.na(res$Nprobes))
	# 		res$Gene.Symbol[idx] <- "---"
	# 		res$Nprobes[idx] <- 1
	# 		res$Possible.Probes[idx] <- res$ID[idx]
	# 	
	# 		n.extra <- nrow(tt.probes)
	# 	}
	# 	else {
	# 		res <- tt.genes
	# 	}
	# }
	# else {
	# 	res <- tt.genes
	# }
	# res <- res[order(as.numeric(res$old_order)), ]
	# rownames(res) <- 1:nrow(res)
	# res$old_order <- NULL
	
	if( verbose ) {
		cat(sprintf("There were %d probes, which collapsed down to %d genes. (%d are unique genes, %d are unannotated probes)\n", n.probes, nrow(res), n.genes, n.extra))
	}
	
	return( res )
}
