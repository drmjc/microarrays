#' Multipanel barplot of results from an lmFit
#' 
#' Very flexible function to barplot results from an lmFit.
#' It handles results from 2 styles of limma analysis:\cr
#' 1: \dQuote{Standard analysis}: model.matrix > lmFit > eBayes > topTable\cr
#' 2: \dQuote{Constrast analysis}: model.matrix > lmFit > fit.constrasts > eBayes >
#' topTable\cr
#' 
#' 1. \dQuote{Standard analysis}\cr
#' for each probe, do a barplot of the normalised data, then
#' an errorbar plot utilising the coefficients
#' and the standard errors (\code{stdev.unscaled * sigma}) from the lmFit1 object.
#' 2. \dQuote{Contrast analysis}\cr
#' for each probe, make 3 barplots.
#' The first 2 are same as standard analysis, the 3rd is an errorbar plot based
#' on fit2 object which you get after doing a \code{\link[limma]{contrasts.fit}}
#' 
#' Improvements\cr
#' Probe selection:\cr
#' 1. probe=a numeric vector of row indices into the lmFit (ie same row order
#' as data)\cr
#' 2. probe=vector of probesetID's which are in the rownames of data and fit1
#' [and fit2] [and calls]\cr
#' 3. supply a topTable object, and set the number of rows from top to bottom
#' to plot.\cr
#' this can be from an F-test or t-test
#' 
#' Colouring array data\cr
#' if you supply a 'calls' object which is same dim as data, and contains \dQuote{P},
#' \dQuote{M} or \dQuote{A}, 
#' then the bars for the expression data will be coloured green, orange or red,
#' respectively.
#' 
#' @param fit1 objects from lmFit. fit2 is optional.
#' @param fit2 objects from lmFit. fit2 is optional. fit1 can also be just a
#'   data.frame which is useful for paired analyses, where you often convert
#'   the expression data (2n) into expression ratios (1n), prior to then
#'   fitting a linear model. In this instance, ags should be fit1=ratios,
#'   fit2=lmFitXYZ, data=rma.
#' @param data data.frame of expression level data
#' @param calls optional data.frame of calls, same dim as rma
#' @param data.type \dQuote{1colour} or \dQuote{2colour}
#' @param probes optional vector of probe indices, or probeset ids
#' @param tt optional toptable of results. if supplied, you should set number
#'   to some
#' @param number optional toptable of results. if supplied, you should set
#'   number to some positive integer corresponding to number of genes to plot.
#' @param probe2genesymbol 2 column table with probe ID's and gene symbols, respectively
#' @param fit1.colour optional vector of colours for the N columns in fit1.
#'   defaults to grey
#' @param fit2.colour optional vector of colours for the N columns in fit1.
#'   defaults to grey
#' @param data.colour optional vector of colours for the N columns in fit1.
#'   defaults to grey. this is ignored if calls != NULL
#' @param hgrid.col do you want horizontal grid lines? NULL means no, otherwise
#'   choose a single colour.
#' @param do.par logical: set the layout and the par settings?
#' @param drop.fit1.intercept logical: drop the intercept term in the first fit object?
#' @param legend.pos Position of the legend. See \code{\link{legend}}, Default = \dQuote{bottomright}
#' @return none
#' @author Mark Cowley, 2009-07-16
#' @examples
#' \dontrun{
#' barplot_lmFit(fit1, data=rma, probes=c("10543233", "10411107"))
#' barplot_lmFit(fit1, fit2, data=rma, probes=c("10543233", "10411107"))
#' barplot_lmFit(fit1, data=rma, tt=topTable, number=2)
#' barplot_lmFit(fit1, fit2, data=rma, tt=topTable, number=2)
#' barplot_lmFit(fit1, fit2, data=rma, calls=calls, tt=topTable, number=2)
#' }
#' @export
barplot_lmFit <- function(
	# lmFit objects
	fit1, fit2=NULL,
	# expression data and P/M/A calls
	data, calls=NULL, 
	data.type=c("1colour", "2colour")[1],
	# probe selection arguments
	probes=NULL, tt=NULL, number=10, 
	probe2genesymbol=NULL,
	# barplot colouring
	fit1.colour="#53406A", fit2.colour="#4F81BD", data.colour=NULL,
	hgrid.col="black", do.par=TRUE,
	# coefficient dropping
	drop.fit1.intercept=FALSE,
	legend.pos="bottomright" ) {

	opar <- par(no.readonly=TRUE)
	if( do.par )
		on.exit(par(opar))

	if( ! data.type %in% c("1colour", "2colour") )
		stop("data.type must be one of \"1colour\" or \"2colour\"")

	if( is.null(probes) && is.null(tt) ) {
		stop("Must supply a vector of probeset id's, or a toptable object.\n")
	}
	else if( is.null(probes) ) {
		probes <- tt$ID[1:number]
	}

	data <- subset(data, rownames(data) %in% probes, drop=FALSE)
	fit1 <- subset.MArrayLM(fit1, probes)
	# if( is.matrix.like(fit1$coef) )
	# 	fit1 <- fit1[probes, 1:ncol(fit1)] # should work even if fit is an lmFit, or a data.frame
	# else
	# 	fit1 <- fit1[probes]

	# try(fit2 <- fit2[probes, ], silent=TRUE)
	if( !is.null(fit2) ) fit2 <- subset.MArrayLM(fit2, probes)
	try(calls <- subset(calls, rownames(calls) %in% probes, drop=FALSE), silent=TRUE)
	try(tt <- tt[match(probes, tt$ID), ], silent=TRUE)
	
	main <- c("normalised data", "model coefficients", "contrast coefficients")
	ylab <- c("Expression level (log2)", "logFC (+/- SE)", "logFC (+/- SE)")
	if( data.type == "2colour" )
		ylab[1] <- "Expression ratio (log2)"

	if( do.par ) {
		mat <- c(ncol(data), ncol(fit1))
		if( !is.null(fit2) ) {
			mat <- c(mat, ncol(fit2))
		}
		layout(matrix(1:length(mat),1), widths=mat+sum(par()$mar[c(2,4)]))
		par(las=2, oma=c(0,0,4,0))
	}	
	
	# nice plot labels. use the top table if supplied
	labels <- paste("ProbeSetID:", probes)
	if( !is.null(tt) ) {
		labels <- topTable.to.label(tt, probes, probe2genesymbol)
	}
	names(labels) <- probes
	
	
	if( !is.null(calls) ) {
		data.colour <- t(apply(calls, 1, calls2colour))
		legend.colours <- calls2colour(c("P", "M", "A"))
	}
	else if( is.null(data.colour) ) {
		data.colour <- as.data.frame(matrix("grey", nrow(data), ncol(data)), stringsAsFactors=FALSE)
		rownames(data.colour) <- probes
		legend.pos <- NA
	}
	else if( is.vector(data.colour) ) {
		data.colour <- as.data.frame(matrix(data.colour, nrow(data), ncol(data), byrow=TRUE), stringsAsFactors=FALSE)
		rownames(data.colour) <- probes
		legend.pos <- NA
	}
	

	if( drop.fit1.intercept && is.matrix(fit1$coefficients) && ncol(fit1$coefficients) > 1 ) {
		rn <- rownames(fit1$coefficients)
		cn <- colnames(fit1$coefficients)
		fit1$coefficients <- fit1$coefficients[,2:ncol(fit1$coefficients)]
		if( is.vector(fit1$coefficients) ) {
			fit1$coefficients <- as.matrix(fit1$coefficients)
			dimnames(fit1$coefficients) <- list(rn, cn[2])
		}
		fit1$stdev.unscaled <- fit1$stdev.unscaled[,2:ncol(fit1$stdev.unscaled)]
		if( is.vector(fit1$stdev.unscaled) ) {
			fit1$stdev.unscaled <- as.matrix(fit1$stdev.unscaled)
			dimnames(fit1$stdev.unscaled) <- list(rn, cn[2])
		}
	}
	
	for(probe in probes) {
		#
		# plot the array data, coloured by calls.
		#
		val <- as.numeric(data[probe,])
		cols <- unlist(data.colour[probe,])
		barplot(val, main=main[1], ylab=ylab[1], names.arg=colnames(data), col=cols)
		if(!is.null(hgrid.col)) {
			hgrid(col=hgrid.col)
		}
		abline(h=0)
		if( !is.na(legend.pos) && legend.pos != "none" ) {
			legend(legend.pos, legend=c("Present", "Marginal", "Absent"), fill=legend.colours, bg="white", inset=0.02, horiz=TRUE)
		}
		
		#
		# plot the first lmFit results.
		#
		if( class(fit1) == "MArrayLM" ) {
			tmp.lmfit <- subset.MArrayLM(fit1, probe)
			# val <- as.numeric(fit1$coefficients[probe,])
			val <- tmp.lmfit$coefficients
			# calc SE see ?lm.series
			# err <- as.numeric(fit1$stdev.unscaled[probe,]) * as.numeric(fit1$sigma[match(probe, rownames(fit1$t))])
			err <- tmp.lmfit$stdev.unscaled * tmp.lmfit$sigma
			# check for NA's in the model fit.
			val[is.na(val)] <- 0
			err[is.na(err)] <- 0
			
			errorbarplot(height=val, errors=err, names.arg=colnames(fit1), 
				main=main[2], ylab=ylab[2], col=fit1.colour, beside=TRUE)#, width=1, space=0.1)
			if(!is.null(hgrid.col)) hgrid(col=hgrid.col)
			abline(h=0)
		}
		else if( is.matrix.like(fit1) ) {
			val <- as.numeric(fit1[probe,])
			val[is.na(val)] <- 0
			barplot(height=val, names.arg=colnames(fit1), 
				main=main[2], ylab=ylab[2], col=fit1.colour)
			if(!is.null(hgrid.col)) hgrid(col=hgrid.col)
			abline(h=0)
			
		}
		else {
			stop("unsupported type for fit1.\n")
		}

		#
		# Optionally plot the 2nd lmFit object.
		#
		if( !is.null(fit2) ) {
			tmp.lmfit <- subset.MArrayLM(fit2, probe)
			
			val <- tmp.lmfit$coefficients
			err <- tmp.lmfit$stdev.unscaled * tmp.lmfit$sigma
			val[is.na(val)] <- 0
			err[is.na(err)] <- 0

			errorbarplot(height=val, errors=err, names.arg=colnames(fit2), 
				main=main[3], ylab=ylab[3], col=fit2.colour, beside=TRUE)#, width=1, space=0.1)
			if(!is.null(hgrid.col)) hgrid(col=hgrid.col)
			abline(h=0)
		}
		mtext(side=3, outer=TRUE, labels[probe], cex=1.5, las=1)
	}
}

