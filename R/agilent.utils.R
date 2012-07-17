#
# Utilities for Agilent arrays
#

###############################################################################
###############################################################################
###############################################################################
# 			IMPORT
###############################################################################
###############################################################################
###############################################################################

#' Import Agilent GeneView files
#' 
#' Import a GeneView file produced by Agilent Feature Extractor.
#' @param file the path to the file to import
#' @author Mark Cowley, 2008-08-01
#' @export
import.Agilent.GeneView.File <- function(file) {
	gv <- read.delim(file, as.is=TRUE, header=TRUE, skip=1)
	rownames(gv) <- gv$SystematicName

	return(gv)
}


#' Import a collection of GeneView files
#' 
#' Import a set of Agilent GeneView files that logically form an experiment.
#' 
#' This first imports each GeneView file, the collates the information into
#' 3 \code{data.frame}s, with ngenes x nsamples, called:\cr
#' \code{value}: the gTotalGeneSignal column\cr
#' \code{error}: the gTotalGeneError column\cr
#' \code{detected}: the gIsDetected column\cr
#' In addition to these, there is a \code{vector} from the \code{ControlType} column (which
#' is identical across all arrays), and a \code{vector} called \code{ProbeIDs}.
#' 
#' @param files the file names. (These will NOT be reordered)
#' @param names if specified, then these names will be used in the column names
#'   of the data.frames; if NULL, then a name based on the filename will be
#'   used
#' @param debug logical
#' @return a list with 5 elements as indicated in Details.
#' @author Mark Cowley, 2008-08-01
#' @export
import.Agilent.GeneView.Experiment <- function(files, names, debug=FALSE) {
	if( missing(names) )
		names <- strip.fileextension( basename(files) )
	else if( length(names) == 1 && names == "auto" ) {
		si <- make.sample.info.Agilent(files, FALSE)
		names <- si$SampleID
	}
	
	raw <- list()
	for(i in 1:length(files)) {
		file <- files[i]
		if( debug ) cat(file, " -> ", names[i], "\n")
		raw[[i]] <- import.Agilent.GeneView.File(file)
	}
	# ensure that the rows are all in the same order...
	probeIDs <- raw[[1]]$SystematicName
	for(i in 2:length(raw)) {
		raw[[i]] <- raw[[i]][match(probeIDs, raw[[i]]$SystematicName), ]
	}
	
	res <- list()	
	res$value <- as.data.frame(lapply(raw, "[", "gTotalGeneSignal"), stringsAsFactors=FALSE)
	res$error <- as.data.frame(lapply(raw, "[", "gTotalGeneError"), stringsAsFactors=FALSE)
	res$detected <- as.data.frame(lapply(raw, "[", "gIsGeneDetected"), stringsAsFactors=FALSE)
	res$controlType <- raw[[1]]$ControlType
	res$probeIDs <- probeIDs
	
	res[1:3] <- lapply(res[1:3], function(x) {rownames(x) <- probeIDs; colnames(x) <- names; x})
	names(res$controlType) <- res$probeIDs

	return( res )
}

# since R does a good job of interpreting what data type each column in a data.frame
# should be, we can safely ignore this...
#
# Mark Cowley, 2008-08-01
#
# .replace.agilent.classes <- function(x) {
# 	x[x=="text"] <- "character"
# 	x[x=="boolean"] <- "logical"
# 	x[x=="float"] <- "double"
# 	# integer -> integer
# 	# 
# }



#' Extract a column from Agilent files.
#' Take a vector of filenames, each of which is a large table from Agilent, 
#' read each table in, then pull out a single column from each table, and
#' stitch these columns together into 1 large table.
#' 
#' @param files a vector of file names, including their paths
#' @param colname the relevant column name to pull out. Beware of colnames with
#'   spaces.
#' @return a table of \code{N} columns, correspnding to the \code{N} files.
#' @author Mark Cowley, 2008-08-29
#' @export
import.agilent.column2table <- function(files, colname="LogRatio") {
	stopifnot( all(file.exists(files)) )
	
	raw <- import.agilent( files[1] )
	colidx <- match(colname, colnames(raw))
	stopifnot( length(colidx) == 1 )

	res <- as.data.frame(matrix(NA, nrow(raw), length(files)))
	
	for( i in 1:length(files) ) {
		tmp <- import.agilent( files[i] )
		res[,i] <- tmp[,colidx]
		assign("tmpres", res, pos=1)
		cat(".")
	}
	cat("\n")
	colnames(res) <- basename(files)
	
	return( res )
}

#' Import an Agilent result file
#' 
#' Import a txt file produced by the Agilent Feature Extraction program.
#' More specifically, import the FEATURE rows, which constitute the 3rd
#' section of the file, which is the majority.
#' 
#' @param f the path to an Agilent Microarray TXT file
#' @author Mark Cowley, 2008-08-29
#' @export
import.agilent <- function(f) {
	if( grepl("gz$", f) ) {
		fout <- tempfile()
		file.gunzip(f, fout)
		on.exit( unlink(fout) )
		f <- fout
	}
	else if( grepl("zip$", f) ) {
		f <- zip.file.extract(f, basename(f))
		# on.exit( unlink(f) )
	}
	
	hdr <- readLines(f, 20)
	skip <- grep("^FEATURES", hdr)[1] - 1
	res <- read.delim(f, skip=skip)[,-1]
	return( res )
}

###############################################################################
###############################################################################
###############################################################################
#			parse the file names to make a sample.info object
###############################################################################
###############################################################################
###############################################################################

#' Autocreate sample.info from Agilent filenames
#' From a vector of Agilent GeneView file names, create a generic \code{sample.info}
#'  object.
#' 
#' @param filenames a character vector of filenames
#' @param outfile if \code{NULL} or \code{FALSE}, then do not write this to a file. if 
#' \code{TRUE}, then write to \dQuote{./sample.info.xls}, otherwise write to the specified file.
#' @return a \code{data.frame} with these columns: \dQuote{FileName}, \dQuote{Barcode}, 
#'   \dQuote{SlideID}, \dQuote{ArrayID}, \dQuote{SampleID}
#' @author Mark Cowley, 2008-08-01
#' @export
make.sample.info.Agilent <- function(filenames, outfile=NULL) {
	stopifnot( all(file.exists(filenames)) )
	
	filenames <- basename(filenames)
	fn <- sub("_GeneView.txt", "", filenames)
	barcode <- sub("_.*", "", fn)
	slideID <- as.numeric(factor(barcode))
	arrayID <- sub("[^_]+_[^_]+_", "", fn)
	sampleID <- paste(slideID, arrayID, sep="_")
	
	sample.info <- data.frame(FileName=filenames, Barcode=barcode, SlideID=slideID, ArrayID=arrayID, SampleID=sampleID, stringsAsFactors=FALSE)
	
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



###############################################################################
###############################################################################
###############################################################################
# 		pre-process and log the data
###############################################################################
###############################################################################
###############################################################################

#' Log transform and offset of Agilent data.
#' 
#' Usually, >90% of the negative values are in [-8,0), with the remainder
#' being < 8.
#' This method adds 8 to all values, then hard thresholds the remaining values
#' that are still negative to some small positive number (eg 1)
#' 
#' @param data a \code{matrix}-like object of agilent data, or a \code{list} 
#' of agilent data with an element called \sQuote{value}
#' @param offset the first offset that is added to all values.
#' @param minVal if any values are small than this minVal, then hard-threshold
#'   them to this minVal
#' @param base the logarithm base. default=2 for log-base-2
#' @return if a matrix-like object was passed into this method, then a log
#'   transformed matrix is returned. If a list was passed into this method,
#'   then the value and error will be logged (base 2).
#' @author Mark Cowley, 2008-08-01
#' @export
transform_Agilent_log <- function(data, offset=8, minVal=1, base=2) {
	if( is.list(data) && "value" %in% names(data) ) {
		data$value <- transform_Agilent_log(data$value, offset=offset, minVal=minVal, base=base)
		data$error <- log(data$error, base=base)
		return( data )
	}
	else {
		res <- data + offset
		negatives <- res < minVal
		if( sum(negatives) > 0 ) {
			cat("Setting", sum(negatives), "values to: ", minVal, "\n")
			res[res < minVal] <- minVal
		} 

		stopifnot( all(res>0) )

		res <- log(res, base=base)
		
		return( res )
	}
}



###############################################################################
###############################################################################
###############################################################################
# 		pre-process and normalise the MJC way
###############################################################################
###############################################################################
###############################################################################

#' Preprocess Agilent miRNA data.
#' From a table of imported Agilent data,
#'  add an offset, truncate, log-transform, normalize and plot a set of 
#' Agilent microarray files. This combines \code{\link{transform_Agilent_log}},
#' \code{\link{agilent.miRNA.normalise.GX10}} or \code{\link[limma]{normalizeBetweenArrays}},
#' and \code{\link{agilent.miRNA.filter.GX10}}
#'
#' @param data a \code{list} of imported Agilent data. see \code{\link{import.agilent}}
#' @param offset add an offset. see \code{transform_Agilent_log}
#' @param min.thresh if any values are < this value, then truncate them to this value. 
#'     see \code{\link{transform_Agilent_log}}
#' @param percentile The percentile to perform percentile-shift normalization. only used if
#'    norm.method=="percentile", see \code{\link{agilent.miRNA.normalise.GX10}}
#' @param min.Pcount integer exclude genes detected in fewer than \code{min.Pcount}
#'    samples. see \code{\link{agilent.miRNA.filter.GX10}}
#' @param species the 3 letter code for which species of miRNA's to include. 
#'    see \code{\link{agilent.miRNA.filter.GX10}}
#' @param plot logical: add boxplots of data during the transformation process
#' @param do.par logical: configure the plotting device settings?
#' @param norm.method The normalisation method to use. one of \dQuote{percentile}, 
#'     \dQuote{none}, \dQuote{scale}, \dQuote{quantile}. If \dQuote{percentile}, then
#'     make sure you set the \code{\link{percentile}} parameter appropriately. If
#' \dQuote{scale} or \dQuote{quantile}, then the \code{\link[limma]{normalizeBetweenArrays}}
#' from \code{limma} is used. If \dQuote{none} then no normalisation is performed.
#' @return a list of log-transformed, filtered, normalised Agilent microRNA data.
#' @export
#' @importFrom limma normalizeBetweenArrays
#' @author Mark Cowley, 2011-08-02
agilent.miRNA.preprocess <- function(data, offset=8, min.thresh=1.0, percentile=0.75, min.Pcount=1, species="mmu", plot=FALSE, do.par=TRUE, norm.method=c("percentile", "none", "scale", "quantile")) {
	if( plot && do.par )
		par(mfrow=c(2,2), mar=c(10,5,4,2)+0.1, las=2)
	
	if( plot ) boxplot(data$value, main="raw levels", ylab="Expression level (unlogged)")
	
	data <- transform_Agilent_log(data, offset, min.thresh)
	if( plot ) boxplot(data$value, main=sprintf("offset by %d, thresholded @ %d and log2 transformed", offset, min.thresh), ylab="Expression level (log2)")

	norm.method <- norm.method[1]
	if( norm.method == "percentile" ) {
		data <- agilent.miRNA.normalise.GX10(data, percentile)
		if( plot ) boxplot(data$value, main=sprintf("%dth percentile shifted", round(percentile*100,0)), ylab="Expression level (log2)")
	}
	else {
		cl <- class(data$value)
		data$value <- normalizeBetweenArrays(as.matrix(data$value), method=norm.method)
		if( cl == "data.frame" )
			data$value <- as.data.frame(data$value, stringsAsFactors=FALSE)
		if( plot ) boxplot(data$value, main=paste("normalised using", norm.method), ylab="Expression level (log2)")
	}

	data <- agilent.miRNA.filter.GX10(data, min.Pcount, species)
	data$Pcount <- rowSums(data$detected)
	
	if( plot ) boxplot(data$value, main=sprintf("filtered on P>=%d and %s miRNA", min.Pcount, species), ylab="Expression level (log2)")

	data
}

###############################################################################
###############################################################################
###############################################################################
# 		pre-process and normalise the GeneSpring GX 10 way
###############################################################################
###############################################################################
###############################################################################


#' Preprocess Agilent miRNA data
#' Preprocess raw Agilent data, then threshold it to a minimum value, log
#' transform, normalise each array such that the 75th percentile is 0.0, throw
#' away miRs from other species, and keep those that are detected at least
#' min.Pcount times. Optionally makes a 2x2 portrait layout plot @@ each step.
#' This is the default procedure by GeneSpring GX 10 (not sure about newer versions)
#' 
#' @param data a \code{list} of imported Agilent data. see \code{\link{import.agilent}}
#' @param min.thresh if any values are < this value, then truncate them to this value. 
#'     see \code{\link{agilent.miRNA.threshold.GX10}}
#' @param percentile The percentile to perform percentile-shift normalization. See
#'    \code{\link{agilent.miRNA.normalise.GX10}}
#' @param min.Pcount integer exclude genes detected in fewer than \code{min.Pcount}
#'    samples. see \code{\link{agilent.miRNA.filter.GX10}}
#' @param species the 3 letter code for which species of miRNA's to include. 
#'    see \code{\link{agilent.miRNA.filter.GX10}}
#' @param plot logical: add boxplots of data during the transformation process
#' @param do.par logical: configure the plotting device settings?
#' @return a set of thresholded, log-transformed, filtered, normalized Agilent
#' miRNA data
#' @author Mark Cowley, 2009-07-22
#' @export
agilent.miRNA.preprocess.GX10 <- function(data, min.thresh=1, percentile=0.75, min.Pcount=1, species="mmu", plot=FALSE, do.par=TRUE) {
	if( plot && do.par )
		par(mfrow=c(2,2), mar=c(10,5,4,2)+0.1, las=2)
	
	if( plot ) boxplot(data$value, main="raw levels", ylab="Expression level (unlogged)")
	
	data <- agilent.miRNA.threshold.GX10(data, min.thresh)
	if( plot ) boxplot(data$value, main=sprintf("thresholded @ %d and log2 transformed", min.thresh), ylab="Expression level (log2)")

	data <- agilent.miRNA.normalise.GX10(data, percentile)
	if( plot ) boxplot(data$value, main=sprintf("%dth percentile shifted", round(percentile*100,0)), ylab="Expression level (log2)")

	data <- agilent.miRNA.filter.GX10(data, min.Pcount, species)
	if( plot ) boxplot(data$value, main=sprintf("filtered on P>=%d and %s miRNA", min.Pcount, species), ylab="Expression level (log2)")

	data
}


#' Threshold Agilent miRNA data
#' Threshold raw Agilent miRNA data such that all values < thresh become =
#' thresh Optionally also log-base-2 transform the results afterwards
#' 
#' @param data a \code{list} of imported Agilent data. see \code{\link{import.agilent}}
#' @param thresh if any values are < this value, then truncate them to this value.
#' @param log logical: perform log-base-2 transformation on the truncated data?
#' @author Mark Cowley, 2009-07-22
#' @export
agilent.miRNA.threshold.GX10 <- function(data, thresh=1.0, log=TRUE) {
	data$value[data$value < thresh] <- thresh
	if( log )
		data$value <- log2(data$value)
	data
}


#' Normalise Agilent miRNA data
#' Normalise an Agilent miRNA data set using the default method from Agilent
#' GeneSpring GX 10 which is percentile shift normalisation. From the GX 10
#' manual:\cr
#' 18.2.1 Percentile Shift Normalization\cr
#' Percentile shift normalization is a global normalization, where the location
#' of all the spot intensities in an array are adjusted. This normalization
#' takes each column in an experiment independently, and computes the
#' percentile of the expression values for this array, across all spots (where
#' n has a range from 0-100 and n=50 is the median). It subtracts this value
#' from the expression value of each entity.
#' 
#' @param data a \code{list} of imported, thresholded, log-base-2 Agilent data.
#'    see \code{\link{agilent.miRNA.threshold.GX10}}
#' @param pctile the percentile to normalize to. default=0.75, ie the 75th 
#'    percentile. Must be in [0,1]
#' @return a list of normalised miRNA data
#' @author Mark Cowley, 2009-07-22
#' @export
agilent.miRNA.normalise.GX10 <- function(data, pctile=0.75) {
	data$value <- colapply(data$value, function(x) {
		x - percentile(x,pctile)
	})
	data$value <- as.data.frame(data$value, stringsAsFactors=FALSE)
	data
}


#' Filter Agilent miRNA data
#' Keep only the miRNA's that start with a certain pattern, such as \dQuote{hsa}, \dQuote{mmu},
#' \dQuote{rno}, etc..., and only those miRNA's that are detected >= min times.
#' @param data a set of normalised Agilent miRNA data. see \code{\link{agilent.miRNA.normalise.GX10}}
#' @param min Retain probes detected in at least \code{min} samples
#' @param species the three letter code for which miRNA's to keep. Set to \code{NULL} 
#'    to include all miRNA's
#' @author Mark Cowley, 2009-07-22
#' @export
agilent.miRNA.filter.GX10 <- function(data, min=1, species="mmu") {
	N <- nrow(data$value)
	if( is.null(species) ) {
		idx.species <- 1:length(data$probeIDs)
	}
	else {
		idx.species <- grep(species, data$probeIDs)
	}
	Pcount <- rowSums(data$detected)
	idx.on <- which(Pcount>=min)
	idx <- intersect(idx.on, idx.species)
	cat(sprintf("Filtering out %d probes. %d were wrong species, and %d were detected < %d times (some could be both).\n", N-length(idx), N-length(idx.species), N-length(idx.on), min))
	data$value <- data$value[idx, ]
	data$error <- data$error[idx, ]
	data$detected <- data$detected[idx, ]
	data$probeIDs <- data$probeIDs[idx]
	data
}


#' cbind Agilent miRNA datasets
#' Combine 2 Agilent miRNA objects side by side. Since Agilent objects are
#' lists of tables, this works on all relevant elements
#' @param x an Agilent miRNA object (which for now is a \code{list})
#' @param y an Agilent miRNA object (which for now is a \code{list})
#' @return a new, larger Agilent miRNA object
#' @export
#' @author Mark Cowley, 2011-08-02
#' @importFrom genomics mirsort
agilent.miRNA.cbind <- function(x, y) {
	probes <- mirsort( union(x$probeIDs, y$probeIDs) )
	.same.order <- function(z, probes) {
		missing.probes <- setdiff(probes, z$probeIDs)
		tmp <- as.data.frame(matrix(NA, length(missing.probes), ncol(z$value)))
		dimnames(tmp) <- list(missing.probes, colnames(z$value))
		z$value <- rbind(tmp,z$value)[probes,]
		z$error <- rbind(tmp,z$error)[probes,]
		z$detected <- rbind(tmp,z$detected)[probes,]
		tmp <- rep(NA, length(missing.probes))
		names(tmp) <- missing.probes
		z$controlType <- c(tmp, z$controlType)[probes]
		z$probeIDs <- probes
		if( "Pcount" %in% names(z) )
			z$Pcount <- c(tmp, z$Pcount)[probes]
		z
	}
	x <- .same.order(x, probes)
	y <- .same.order(y, probes)
	
	res <- list()
	res$value <- cbind(x$value, y$value)
	res$error <- cbind(x$error, y$error)
	res$detected <- cbind(x$detected, y$detected)
	res$controlType <- x$controlType
	res$controlType[is.na(x$controlType)] <- y$controlType[is.na(x$controlType)]
	res$probeIDs <- probes
	if( ("Pcount" %in% names(x)) && ("Pcount" %in% names(y)) ) {
		res$Pcount <- x$Pcount
		res$Pcount[is.na(x$Pcount)] <- 0
		tmp <- y$Pcount
		tmp[is.na(y$Pcount)] <- 0
		res$Pcount <- res$Pcount + tmp
	}
	
	res
}
