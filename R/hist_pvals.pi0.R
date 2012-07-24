#' histrogram of p-values
#' 
#' Plot a histogram of unadjusted P-values, and the expected proportion of
#' truly nulls, and the estimated pi0 from fitting the non-\code{NA} pvalues through
#' qvalue with default options.
#' 
#' @param pvals a numeric vector of unadjusted pvalues, in [0,1]
#' @param bins the number of histogram bins
#' @param main the plot title
#' @param xlab the plot x-axis label
#' @param annotate logical: if \code{TRUE}, annotate the plots with 2 horizontal
#' lines, at pi0 under H0, and estimated pi0
#' @param theme one of \dQuote{red} or \dQuote{black}, where \dQuote{red} colours the bars
#' from deep red to white, or \dQuote{black} leaves the bars black & white.
#' 
#' @return nothing.
#' 
#' @author Mark Cowley, 2008-10-02
#' @export
#' @importFrom mjcstats qvalue2
#' 
hist_pvals_pi0 <- function(pvals, bins=20, main="", xlab="unadjusted P-values", annotate=TRUE, theme=c("red", "black")) {
	theme <- theme[1]
	breaks <- seq(0,1,length.out=bins+1)
	
	if(theme == "red") {
		col <- colours.mjc("reds", bins)
		border <- NA
	}
	else if(theme == "black") {
		col <- NULL
		border <- NULL
	}
		
	hist(pvals, breaks=breaks, main=main, xlab=xlab, col=col, border=border)
	if( annotate ) {
		q <- qvalue2(na.rm(pvals)) # @TODO, in case of lots of NA's, shouldn't this force n=length(pvals)
		abline(h=length(pvals)/bins, lty="dashed")
		abline(h=length(pvals)/bins*q$pi0, lty="dotted")
		legend("topright", c("pi0 under H0", "estimated pi0"), lty=c("dashed", "dotted"), inset=0.02)
	}
}
