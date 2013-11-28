#' plot_exprs_vs_detection
#' 
#' From a set of expression data, make one 2-panel-plot for each probe in \code{probes},
#' where the top panel is the sorted expression levels, from low to high, and the bottom
#' panel is the detection p-values. Datapoints are \dQuote{x} or \dQuote{o} depdending on
#' whether they are above or below the \code{detection.thresh}. Optionally, certain named
#' \code{samples} can the highlighted in red. A percentile axis is added to the top of
#' the top panel. If you supply more than 1 probe, then multiple plots will be displayed.
#' If you leave \code{mains=missing}, then you must specify \code{gene.column}, as the column
#' name within \code{fData(x)[,gene.column]} which contains gene symbols.
#' It's been designed to display large numbers of samples.
#'
#' @param x an ExpressionSet object
#' @param probes a character vector of at least 1 probe ID. must be in \code{featureNames(x)}
#' @param samples an optional vector of sample names to highlight in the plot
#' @param detection.thresh an optional detection threshold for calling 'detected' genes
#' @param mains an optional plot title(s). If missing, then the plot title is set to
#'  gene - probe - n=<n>, where gene is from the \code{fData(x)[,gene.column]}, probe
#'  is from \code{probes}, and <n> is the number of samples in \code{x}.
#' @param gene.column the column name within \code{fData(x)} which contains the gene symbols.
#' 
#' @return nothing; a 2-panel plot is produced
#' 
#' @author Mark Cowley
#' @exportMethod plot_exprs_vs_detection
#' @rdname plot_exprs_vs_detection-methods
#' @docType methods
#' @examples
#' if( require(lumi) ) {
#'   data(example.lumi)
#'   example.lumi
#'   plot_exprs_vs_detection(example.lumi, probes="oZsQEQXp9ccVIlwoQo", samples=c("A01", "A02"), 0.01, gene.column="TargetID")
#' }
setGeneric(
	"plot_exprs_vs_detection",
	function(x, probes, samples, detection.thresh=0.01, mains="Expression vs Detection", gene.column="SymbolReannotated") {
		standardGeneric("plot_exprs_vs_detection")
	}
)

#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,character,character,numeric,character,NULL-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "character", "character", "numeric", "character", "NULL"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		plot_exprs_vs_detection(x, probes, samples, detection.thresh, mains, "SymbolReannotated")
	}
)

#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,character,character,numeric,character,character-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "character", "character", "numeric", "character", "character"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		stopifnot(length(probes) > 0, all(probes %in% featureNames(x)))
		stopifnot(length(samples) == 0 || all(samples %in% sampleNames(x)))
		length(mains) == 0 && length(gene.column) == 0 && stop("Must supply either 'mains' or a 'gene.column'")
		
		x <- x[probes, ]
		
		if( length(mains) == 0 ) { 
			genes <- fData(x)[, gene.column[1]]
			mains <- sprintf("%s - %s - n=%d", genes, probes, ncol(x))
		}
		else {
			mains <- recycle(mains, length(probes))
		}
		
		PMIN <- 0.0001 # minimum detection p-value

		opar <- par(no.readonly=TRUE)
		on.exit(par(opar))

		par(mfrow=c(2,1), mgp=c(3.1,1,0)) 
		for(i in seq(along=probes)) {
			probe <- probes[i]
			main <- mains[i]

			x.subset <- x[probe,]
			x.subset <- x.subset[,order(exprs(x.subset))]

			pch <- c(4, 1)[as.numeric(detection(x.subset) < detection.thresh) + 1]
			col <- rep("black", length(sampleNames(x.subset)))

			if( length(samples) > 0 ) {
				col[match(samples, sampleNames(x.subset))] <- "red"
			}

			par(mar=c(0.5,4.1,4.1,1.1))
			plot(exprs(x.subset)[1,], xlab="", ylab="Expression Level (log2)", main="", pch=pch, col=col, xaxt="n")
			title(main=main, line=3)
			axis.percentiles(side=3, max=length(sampleNames(x.subset)))
			grid.percentiles(side=3, max=length(sampleNames(x.subset)))
			if( length(samples) > 0 ) {
				idx <- match(samples, sampleNames(x.subset))
				# points(idx, exprs(x.subset)[1,idx], pch=19, col="red")
				text(idx, exprs(x.subset)[1,idx] + 0.1*dRange(exprs(x.subset), na.rm=TRUE), samples, col="red", pos=2)
			}

			#
			# detection p-values
			# 
			p <- pmax(detection(x.subset)[1,], PMIN)
			par(mar=c(4.1,4.1,0.5,1.1))
			plot(p, xlab="", ylab="detection pval", log="y", ylim=c(min(c(p, 0.05), na.rm=TRUE), 1), pch=pch, col=col)
			# if( length(samples) > 0 ) {
			# 	idx <- match(samples, sampleNames(x.subset))
			# 	points(idx, p[idx], pch=19, col="red")
			# }
			abline(h=detection.thresh, lty="dashed")
		}
	}
)

#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,missing,missing,missing,missing,missing-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "missing", "missing", "missing", "missing", "missing"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		plot_exprs_vs_detection(x, probes=featureNames(x)[1], samples=character(), detection.thresh=0.01, mains="Expression vs Detection", gene.column=NULL)
	}
)

#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,character,missing,missing,missing,missing-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "character", "missing", "missing", "missing", "missing"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		plot_exprs_vs_detection(x, probes, samples=character(), detection.thresh=0.01, mains="Expression vs Detection", gene.column=NULL)
	}
)


#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,character,character,missing,missing,missing-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "character", "character", "missing", "missing", "missing"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		plot_exprs_vs_detection(x, probes, samples, detection.thresh=0.01, mains="Expression vs Detection", gene.column=NULL)
	}
)


#' @rdname plot_exprs_vs_detection-methods
#' @aliases plot_exprs_vs_detection,ExpressionSet,character,character,numeric,missing,missing-method
setMethod(
	"plot_exprs_vs_detection",
	signature=signature("ExpressionSet", "character", "character", "numeric", "missing", "missing"),
	function(x, probes, samples, detection.thresh, mains, gene.column) {
		plot_exprs_vs_detection(x, probes, samples, detection.thresh, mains="Expression vs Detection", gene.column=NULL)
	}
)

