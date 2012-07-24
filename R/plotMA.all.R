#' Make N MAplots, writing out files into a directory.
#' 
#' @param MA an MAList object
#' @param dir the output directory
#' @param width arguments to png
#' @param height arguments to png
#' @param \dots further args passed to plotMA.
#' @return none. creates N files, one per array.  See also: limma:plotMA
#' @author Mark Cowley, 2009-11-02
#' 
#' @export
#' @importFrom limma plotMA
#' 
plotMA.all <- function(MA, dir, width=1024, height=768, ...) {
	for(i in 1:ncol(MA)) {
		f <- file.path(dir, paste(colnames(MA)[i], ".png", sep=""))
		png(f, width, height)
		plotMA(MA, i, ...)
		dev.off()
	}
}
