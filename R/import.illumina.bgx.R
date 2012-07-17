#' Import an Illumina microarray manifest file.
#' 
#' Illumina microarray manifest files describe the contents of each microarray, 
#' including the probe names and sequences among many other things. You can
#' download them from [1] or [2]. You should download the text versions of 
#' these files, unzip them and import using this function.
#' You can either import the probe information,or the control information.
#' 
#' @note This seems like a way to unzip BGX manifest files which i don't remember
#' being possible:\cr
#'  if the file ends in zip, then unzip it.
#' then unzip the bgx file using \code{gunzip -S .bgx xxxxx}, then rename it to 
#' end with \code{.bgx} again & import it using this function. 
#' 
#' @param file the path to the Illimina BGX file
#' @param what what probes to import? one of \dQuote{probes} or \dQuote{controls}
#' @author Mark Cowley, 2008-10-23
#' @return a \code{data.frame} representation of the probes from a BGX file.
#' @references 
#' [1] \url{http://www.switchtoi.com/annotationfiles.ilmn}
#' [2] \url{http://www.switchtoi.com/annotationprevfiles.ilmn}
#' @export
import.illumina.bgx <- function(file, what=c("probes", "controls")) {
	what <- what[1]
	
	header <- readLines(file, 20)
	nProbes <- as.numeric( sub("^.*\t", "", grep("Number of Probes", header, value=TRUE)) )
	nControls <- as.numeric( sub("^.*\t", "", grep("Number of Controls", header, value=TRUE)) )
	
	headerLength <- which(header == "[Probes]") - 1
	probesFirstRow <- headerLength + 2
	controlsFirstRow <- probesFirstRow + 1 + nProbes + 1
	
	if( what == "probes" ) {
		res <- read.delim(file, skip=probesFirstRow-1, nrows=nProbes)
	}
	else {
		res <- read.delim(file, skip=controlsFirstRow-1, nrows=nControls)
	}
	
	return( res )
}
