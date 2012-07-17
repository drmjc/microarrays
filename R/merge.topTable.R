#' Function to merge N topTable's together, retaining a particular column from
#' each topTable.
#' 
#' @param tt.list a named list of topTable objects.
#' @param keep the column to keep or retain. It should be one of the column
#'   names in each GSEA object.
#' @return a data.frame of N+1 columns: ID, then the values that were extracted
#'   from each topTable. eg the t-stats
#' @author Mark Cowley, 2009-12-16
#' @export
merge_topTable <- function(tt.list, keep=c("t", "adj.P.Val", "P.Value")) {
	keep <- keep[1]
	stopifnot(keep %in% colnames(tt.list[[1]]))
	
	tmp <- lapply(tt.list, function(x) x[,1])
	genes <- sort(unique(unlist(tmp)))

	res <- as.data.frame(matrix(NA, nrow=length(genes), ncol=length(tt.list)+1), stringsAsFactors=FALSE)
	colnames(res) <- c("ID", names(tt.list))
	res$ID <- genes

	for(i in 1:length(tt.list)) {
		idx <- match(tt.list[[i]]$ID, genes)
		res[idx,i+1] <- tt.list[[i]][,match(keep, colnames(tt.list[[i]]))]
	}

	res
}
