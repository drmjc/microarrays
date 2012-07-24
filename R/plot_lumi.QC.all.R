#' Perform all of the QC plots that are mentioned in the lumi vignette.
#' It creates a density plot, boxplot, PCA, HCL, CV-plot, pairs and MA plots.
#' The pairs and MA plot can take a LONG time (>10mins?) if you have lots of arrays,
#' so you can skip these plots using the \code{MA=FALSE} and \code{pairs=FALSE} 
#' arguments
#' 
#' @note The CV plot some times fails, even if you don't have missing data in your LumiBatch. In this case,
#' you'll get this Error:\cr
#' \code{Error in density.default(newX[, i], ...) : 'x' contains missing values}
#' 
#' @param x a LumiBatch object.
#' @param dir the directory to create the plot files
#' @param prefix prefix the plot filename prefix. eg \dQuote{unnorm}, \dQuote{qnorm}, \dQuote{vst-transformed}
#' @param title.prefix default=\dQuote{Unnormalised}
#' @param MA logical: create MA plots vs the average array?
#' @param pairs logical: createa pairs plot? takes ages if you have lots of arrays.
#' 
#' @author Mark Cowley, 2008-10-23
#' @export
#' @importFrom lumi plotSampleRelation MAplot pairs plot
#' @importClassesFrom lumi LumiBatch
#' 
#' @examples
#' \dontrun{
#' dir.create("QC/01.unnorm")
#' plot_lumi_QC_all(x.raw, "QC/01.unnorm/", "raw", "Unnormalised")
#' dir.create("QC/02.transformed")
#' plot_lumi_QC_all(x.transformed, "QC/02.transformed/", "vst", "VST
#' Transformed")
#' dir.create("QC/03.rsn")
#' plot_lumi_QC_all(x.norm.rsn, "QC/03.rsn/", "rsn", "RSN Normalised")
#' }
plot_lumi_QC_all <- function(x, dir, prefix, title.prefix="Unnormalised", MA=TRUE, pairs=TRUE) {
    !missing(x) && is(x, "LumiBatch") || stop("x must be a LumiBatch object")
    !missing(dir) || stop("Must specify dir")
    !missing(prefix) || stop("Must specify prefix")
    !missing(title.prefix) || stop("Must specify title.prefix")
	
    dir.create(dir, showWarnings=FALSE)
	
	message("Density plot")
	N <- dim(x)[2]
	f <- file.path(dir, paste0(prefix, ".density.png"))
	png.SVGA(f)
	par(las=1)
	density(x, main=paste(title.prefix, "density plot"))
	dev.off()
	
	message("Boxplot")
	f <- file.path(dir, paste0(prefix, ".boxplot.png"))
	png.SVGA(f)
	par(las=1)
	boxplot(x, main=paste(title.prefix, "boxplot"))
	dev.off()
	
	minsz <- max(1600, (N+1)*150) # the pairs plots have N+1 x N+1 plots in one.
	if( pairs ) {
		message("xy-pairs plot")
		f <- file.path(dir, paste0(prefix, ".pairs.png"))
		png(f, minsz, minsz)
		par(las=1)
		pairs(x, main=paste(title.prefix, "pairs plot"))
		dev.off()
	}
	
	if( MA ) {
		message("MA plot")
		f <- file.path(dir, paste0(prefix, ".MAplot.png"))
		png(f, minsz, minsz)
		par(las=1)
		MAplot(x, main=paste(title.prefix, "MA plot"))
		dev.off()
	}

	try({
		message("CV plot")
		f <- file.path(dir, paste0(prefix, ".CV.png"))
		png.SVGA(f)
		par(las=1)
		plot(x, what="cv")
		dev.off()
	}, silent=TRUE)

	message("PCA plot")
	f <- file.path(dir, paste0(prefix, ".PCA.png"))
	png.SVGA(f)
	par(las=1)
	plotSampleRelation(x, method="mds", color=rep(c("blue", "red"), each=3))
	dev.off()

	message("HCL plot")
	f <- file.path(dir, paste0(prefix, ".HCL.png"))
	png.SVGA(f)
	par(las=1)
	plot(x, what="sampleRelation")
	dev.off()
	
	message("rank vs stdev plot")
	f <- file.path(dir, paste0(prefix, ".rank-vs-stdev.png"))
	png.SVGA(f)
	par(las=1)
	plot.rank.vs.sd(x, main=paste(title.prefix, "rank vs stdev"))
	dev.off()
	
}
