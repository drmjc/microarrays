#' Sumarise an Affymetrix DABG file
#' The Detected Above Background method by Affymetrix generates a P-Value for the probability
#' of each probeset being expressed above the matched control probes.
#' Import an Affymetrix DABG results table, and: boxplot of these p-values for each array
#' summary table of how many probesets pass certain P-value thresholds. see
#' also \code{\link{plot_dabg_vs_rma}}
#' @param file the path to an Affymetrix DABG result file.
#' @author Mark Cowley, 17/1/08
#' @examples
#' \dontrun{
#' dabg <- summarise.dabg("./dabg.summary.txt")
#' }
#' @export
summarise.dabg <- function(file) {
	dabg <- read.delim(file, as.is=TRUE, comment.char="#", row.names=1)
	colnames(dabg) <- sub(".CEL", "", colnames(dabg))
	# DABG BOXPLOTs
	png(file.path(dirname(file), "dabg.summary.png"), 1400, 800)
	par(mfrow=c(1,2))
	boxplot(dabg, main="Boxplots of DABG", ylab="P-value")
	boxplot(dabg + 1e-06, log="y", main="Boxplots of DABG", ylab="P-value (log scale)")
	dev.off()
	# how many pass each threshold?
	
	
	thresholds <- c(0.05, 0.01, 0.001, 0.0001, 0.00001)
	res <- matrix(NA, length(thresholds)+1, ncol(dabg))
	dimnames(res) <- list(c("P > 0.05", paste("P <", thresholds)), colnames(dabg))

	for(t in 1:length(thresholds)) {
		for(array in 1:ncol(dabg)){
			if( t == 1 ) # also determine how many have P > 0.05
				res[t,array] <- sum(dabg[,array] >= 0.05)
			res[t+1,array] <- sum(dabg[,array] < thresholds[t])
		}
	}

	res <- res / nrow(dabg) * 100
	res <- cbind(res, average=rowMeans(res))
	res <- colapply(res, prettyNum)
	write.delim(res, file.path(dirname(file), "dabg.summary.tsv"), row.names="threshold")
	
	invisible(dabg)
}


#' Density plot of DABG vs expression levels
#' Plot the RMA expression levels as density plots, split by various DABG
#' P-value thresholds.
#' 
#' @param dabg a \code{data.frame} of DABG P-Values
#' @param rma a \code{data.frame} of expression levels (usually RMA normalised levels, thus the name)
#' @param thresholds a vector of DABG p-value thresholds
#' @author Mark Cowley, 17/1/08
#' @examples
#' \dontrun{
#' plot_dabg_vs_rma(dabg, rma)
#' }
#' @export
#' @importFrom stats density
plot_dabg_vs_rma <- function(dabg, rma, thresholds=c(0.05, 0.01, 0.001, 0.0001, 0.00001) ) {
	dabg <- dabg[intersect(rownames(dabg), rownames(rma)),]
	rma <- rma[intersect(rownames(dabg), rownames(rma)),]
	
	auto.mfrow(ncol(dabg))
	for(array in 1:ncol(dabg)) {
		plot(density(rma[,array]), xlim=c(0,15), ylim=c(0, 0.30), main=colnames(dabg)[array], xlab="RMA expression value", col="black")
		lines(density(rma[dabg[,array] >= 0.05, array]), col="darkgrey")
		for(i in 1:length(thresholds))
			lines(density(rma[dabg[,array] < thresholds[i],array]), col=i+1)
	
		legend("topright", c("all", "P >= 0.05", paste("P <", thresholds)), col=c("black", "darkgrey", 1:length(thresholds)+1), bty="n", inset=0.01, lty=1)
	}
	auto.mfrow(ncol(dabg), setup=FALSE)
}
