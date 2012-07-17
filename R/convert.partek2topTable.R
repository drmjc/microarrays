#' Convert Partek to topTable
#' Convert a toptable made by Partek's ANOVA into a limma-like-topTable.
#'  I expect there
#' to be a few columns of annotation, including the Probeset.ID; then there
#' should be at least a p value, t stat, fold change column. If there's a q
#' value column then great it will be used, otherwise q-values will be
#' generated for you. Do NOT let there be > 1 pval/qval/t-stat column...
#' 
#' @param partek a \code{data.frame} created by Partek, and imported into \R
#' @return a \code{data.frame}
#' @author Mark Cowley, 2009-12-21
#' @export
#' @importFrom mjcstats p.adjust
convert.partek2topTable <- function(partek) {
	ID.col <- grep("Probeset.ID", colnames(partek)); stopifnot(length(ID.col) == 1)
	t.col <- grep("^T\\.", colnames(partek)); stopifnot(length(t.col) == 1)
	p.col <- grep("^p\\.value\\.", colnames(partek)); ; stopifnot(length(p.col) == 1)
	fdr.col <- grep("qvalue", colnames(partek))
	if( length(fdr.col) == 0 ) {
		partek$Q <- p.adjust(partek[,p.col], method="q")
		fdr.col <- grep("^Q$", colnames(partek))
	}
	FC.col <- grep("^Fold\\.Change\\.", colnames(partek))[1]; stopifnot(length(FC.col) == 1)
	
	tt <- partek[,c(ID.col, FC.col, t.col, p.col, fdr.col)]
	colnames(tt) <- c("ID", "FC", "t", "P.Value", "adj.P.Val")
	tt$logFC <- log2(abs(tt$FC))*sign(tt$FC)
	tt$B <- NA
	tt <- tt[, c("ID", "logFC", "t", "P.Value", "adj.P.Val", "B")]
	tt <- tt[order(abs(tt$t), decreasing=TRUE), ]
	tt
}
