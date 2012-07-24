# From the results of a RankProd analysis, create a signed score suitable for
# something like GSEA analysis.
#
# Parameters:
#	rp: a rankprod result, or a list containing a table called \dQuote{RPs}
#	method: one of \dQuote{inv.log}, or \dQuote{inverse}
#
# Value:
# a named vector of rank product signed scores
#
# Mark Cowley, 2009-01-09


#' From the results of a RankProd analysis, create a signed score suitable for
#' something like GSEA analysis.
#' 
#' @param rp a rankprod result, or a list containing a table called
#'   \dQuote{RPs}
#' @param method one of \dQuote{inv.log}, or \dQuote{inverse}
#' @return a named vector of rank product signed scores
#' @author Mark Cowley, 2009-01-09
#' @export
RankProd.signed.score <- function(rp, method=c("inv.log", "inverse")) {
	method <- method[1]

	if( is.list(rp) && "RPs" %in% names(rp) )
		rp <- rp$RPs

	which.direction <- apply(rp, 1, which.min)

	if( method == "inverse")
		score <- ifelse(which.direction==1, (1/(rp[,1]))*-1, (1/(rp[,2]))*1) * 1000
	else if( method == "inv.log" )
		score <- ifelse(which.direction==1, (1/log(rp[,1]))*-1, (1/log(rp[,2]))*1) * 1
	else
		stop("Unsupported method.\n")

	if( !is.null(rownames(rp)) )
		names(score) <- rownames(rp)

	score
}


#' CAT plot comparing N RankProd analyses.
#' 
#' All pair-wise combinations of CAT plots are made between each of the
#' supplied lists produced by RP or RPAdvance.
#' 
#' @param \dots at least 1 rankprod object
#' @param names the names of each rankprod object. length >= number of \dots
#'   arguments
#' @return none. a CAT plot is made
#' @author Mark Cowley, 2008-12-10
#' @export
catplot.RankProd <- function(..., names=LETTERS) {
	par(mfrow=c(1,2))

	RP <- list(...)
	N <- length(RP)
	
	col <- 0
	legend.text <- NULL
	for(i in 1:(N-1)) {
		for(j in (i+1):N) {
			col <- col + 1
			new.plot <- FALSE
			if( i == 1 && j == 2 )
				new.plot <- TRUE

			for(direction in c(1,2)) {
				par(mfg=c(1,direction))
				idxA <- order(RP[[i]]$RPs[, direction], decreasing=FALSE)
				idxB <- order(RP[[j]]$RPs[, direction], decreasing=FALSE)
				catplot( rownames(RP[[i]]$RPs)[idxA],
						  rownames(RP[[j]]$RPs)[idxB],
						  col=col, add=!new.plot, 
						  main=c("UP-regulated", "DOWN-regulated")[direction] )
			}
			legend.text <- c(legend.text, paste(names[i], "vs", names[j]))
		}
	}
	legend("bottomright", col=1:col, pch=1, legend.text, inset=0.01)
}
