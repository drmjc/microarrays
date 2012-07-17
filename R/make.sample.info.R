# From a vector of CEL file names, create a generic sample.info object.
#
# Parameters:
#	files: a character vector of filenames
#	outfile: if NULL or FALSE, then do write this to a file
#			 if TRUE, then write to ./sample.info.xls, otherwise write to the specified file.
#	auto.batch: try to determine any processing batches. see batch.thresh
#	batch.thresh: if auto.batch, then increment the numeric batch ID if the time between
#			two consecutive files is greater than this threshold in MINTUES. Default is 1 day/
#
# Value:
#	a data.frame
#
# Mark Cowley, 2008-07-29


#' From a vector of CEL file names, create a generic sample.info object.
#' 
#' @param files a character vector of filenames
#' @param outfile if NULL or FALSE, then do write this to a file if TRUE, then
#'   write to ./sample.info.xls, otherwise write to the specified file.
#' @param auto.batch try to determine any processing batches. see batch.thresh
#' @param batch.thresh if auto.batch, then increment the numeric batch ID if
#'   the time between two consecutive files is greater than this threshold in
#'   MINTUES. Default is 1 day/
#' @return a data.frame
#' @author Mark Cowley, 2008-07-29
#' @export
make.sample.info.Affymetrix <- function(files=NULL, outfile=NULL, auto.batch=TRUE, batch.thresh=1440) {
	if( is.null(files) && file.exists("CEL")) {
		files <- dir("CEL", pattern="CEL$", full.names=TRUE)
	}
	stopifnot( all(file.exists(files)) )
	ID <- strip.fileextension(basename(files))
	dates <- celDate(files)
	sample.info <- data.frame(CEL=basename(files), ID=ID, Param1=NA, Date=dates, stringsAsFactors=FALSE)

	if( auto.batch ) {
		# convert the dates to a format that R understands, then find out the
		# number of mins between each CEL creation date, and some time in 1990,
		# and then the number of mins between each CEL creation dates.
		time0 <- strptime("01/01/90 01:01:01", format="%m/%d/%y %H:%M:%S")
		times <- strptime(dates, format="%m/%d/%y %H:%M:%S")
		dTime <- as.numeric(difftime(times, time0, units="mins"))
		o <- order(dTime, decreasing=FALSE)
		sample.info <- sample.info[o,]
		dTime <- dTime[o]
		times <- times[o]
		ddTime <- c(0, dTime[2:length(dTime)] - dTime[1:(length(dTime)-1)])
		# make a numeric batch ID, which increments if the period between 
		# 2 arrays being processed is greater than the threshold.
		newBatch <- ddTime > batch.thresh
		batch <- cumsum(newBatch) + 1
		sample.info$ScanOrder <- 1:nrow(sample.info)
		sample.info$Batch <- batch
		sample.info <- sample.info[order(sample.info$CEL), ]
	}

	if( !is.null(outfile) && outfile != FALSE ) {
		if( outfile == TRUE ) {
			outfile <- "sample.info.xls"
		}
		if( file.exists(outfile) ) {
			stop("outfile", sQuote(outfile), "exists, please remove it, and re-run this function\n")
		}
		write.delim(sample.info, outfile, na="")
	}


	sample.info
}
