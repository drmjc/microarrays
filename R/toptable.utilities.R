#' Summarise a topTable object
#' 
#' Get the DE genes from a \code{topTable} or \code{topTableQ} object,
#' looking in the \dQuote{P.Value}, and \dQuote{q} or \dQuote{logFC} column,
#' returning
#' a \code{data.frame} with 1 row, with numbers of DE genes in each column.
#' 
#' @param tt a \code{\link[limma]{topTable}} or \code{topTableQ} object
#' @param p.thresh a numberic vector of p-value threhshlds to use. set to NULL
#'   if not required.
#' @param q.thresh a numeric vector of q-value threhsolds to use. set to NULL
#'   if not required
#' @param logFC.thresh a numeric vector of absolute logFC threhsolds to use.
#'   set to NULL if not required
#' @return a 1 x n \code{data.frame}, where n = np + nq where np is the number of
#'   p-value threhsolds, and nq is the number of q-value thresholds.  See also:
#'   \code{\link{export.DEgenes.topTable}}
#' @author Mark Cowley, 25/2/2008
#' @export
summarise.topTable <- function(tt, p.thresh=c(0.05, 0.001, 0.0001), q.thresh=c(0.25, 0.10, 0.05), logFC.thresh=c(0.585,1,2)) {

	FC.thresh <- round(2^logFC.thresh,1)
	
	if( is.topTable.list(tt) ) {
		if( is.null(names(tt)) ) names(tt) <- paste("tt", 1:length(tt), sep="")

		tmp <- lapply(tt, summarise.topTable, 
			p.thresh=p.thresh, q.thresh=q.thresh, logFC.thresh=logFC.thresh)
		res <- rbind.list(tmp)
		res <- data.frame(topTable=names(tt), res, check.names=FALSE)
		rownames(res) <- NULL
		return( res )
	}
	else {
		HAS.LOGFC <- "logFC" %in% colnames(tt)

		# P value summary
		res <- as.data.frame(matrix(NA, 1,1+length(p.thresh)+length(q.thresh)+length(logFC.thresh)), stringsAsFactors=FALSE)
		colnames(res) <- c( paste0("P<",p.thresh), paste0("q<", q.thresh), paste0("FC>", FC.thresh), "N" )
		
		res$N <- nrow(tt)
		
		for(j in 1:length(p.thresh)) {
			res[1,j] <- sum(tt[,"P.Value"] < p.thresh[j], na.rm=TRUE)
		}

		# FDR summary
		qcol <- which(colnames(tt) %in% c("q", "q (Storey)"))
		if( length(qcol) == 0 )
			qcol <- which(colnames(tt) == "adj.P.Val")
		for(j in 1:length(q.thresh)) {
			res[1,j+length(p.thresh)] <- sum(tt[,qcol] < q.thresh[j], na.rm=TRUE)
		}
		
		# logFC summary
		if( HAS.LOGFC ) {
			for(j in 1:length(logFC.thresh)) {
				res[1,j+length(p.thresh)+length(q.thresh)] <- sum(abs(tt[,"logFC"]) > logFC.thresh[j], na.rm=TRUE)
			}
		}
		
		res <- res[,c(ncol(res), 1:(ncol(res)-1))]
		res
	}
}

#' Summarise a list of topTables.
#' 
#' This represnts a new style of summary, where N top tables makes N columns
#' (plus some descriptive ones...)
#' 
#' @param tt.list must be a list of toptable objects, and no F stat tables. all
#'   need to have logFC columns.
#' @param p.thresh a numberic vector of p-value threhshlds to use. set to NULL
#'   if not required.
#' @param q.thresh a numeric vector of q-value threhsolds to use. set to NULL
#'   if not required
#' @param logFC.thresh a numeric vector of absolute logFC threhsolds to use.
#'   set to NULL if not required
#' @return a \code{data.frame} See also: \code{\link{summarise.topTable}},
#'   \code{\link{summarise.topTable.updown}}
#' @author Mark Cowley, 2009-12-21
#' @export
summarise.topTable.list <- function(tt.list, p.thresh=c(0.05, 0.001, 0.0001), q.thresh=c(0.25, 0.10, 0.05), logFC.thresh=c(0.585,1,2)) {

	FC.thresh <- round(2^logFC.thresh,1)
	
	if( !is.topTable.list(tt.list) ) {
		stop("Must be a list of topTable objects.\n")
	}

	tmp <- lapply(tt.list, summarise.topTable.updown, p.thresh=p.thresh, q.thresh=q.thresh, logFC.thresh=c(0,1))
	res <- as.data.frame(matrix(NA, nrow(tmp[[1]]$up)*2, length(tt.list)+2), stringsAsFactors=FALSE)
	colnames(res) <- c("threshold", "direction", names(tt.list))
	res$threshold <- rep(rownames(tmp[[1]]$up), each=2)
	res$direction <- c("up", "down")
	for(i in 1:length(tt.list)) {
		res[seq(1,nrow(res), 2),i+2] <- tmp[[i]]$up[,1]
		res[seq(2,nrow(res), 2),i+2] <- tmp[[i]]$down[,1]
	}
	res
}


#' Summarise a toptable into the up/down regulated genes
#'
#' @param tt a limma \code{topTable}, or \code{list(topTable)}
#' @param p.thresh a vector of p value thresholds, set to \code{NULL} to ignore
#' @param q.thresh a vector of adj.P.Val thresholds, set to \code{NULL} to ignore.
#' @param logFC.thresh a vector of absolute logFC value thresholds
#' @return a data.frame of counts of differentially expressed genes
#' @export
#' @author Mark Cowley
#' @seealso \code{\link{summarise.topTable}}, \code{\link{summarise.topTable.list}}
summarise.topTable.updown <- function(tt, p.thresh=c(0.05, 0.001, 0.0001), q.thresh=c(0.25, 0.10, 0.05), logFC.thresh=c(0.585,1,2)) {
	
	.summ <- function(stats, stat.thresh, lfc, lfc.thresh, dir=c("up", "down", "either") ) {
		dir <- dir[1]
		res <- matrix(NA, length(stat.thresh), length(lfc.thresh))

		# apply some logic to the lFC to avoid 3 if statements in the forfor loop
		if( dir == "down" )
			lfc <- lfc * -1
		else if( dir == "either" )
			lfc <- abs(lfc)

		for( i in seq(along=stat.thresh) ) {
			for( j in seq(along=lfc.thresh) ) {
				res[i,j] <- sum(stats < stat.thresh[i] & lfc > lfc.thresh[j], na.rm=TRUE)
			}
		}
		res
	}
	
	res <- list()
	for(dir in c("up", "down", "either")) {
		p <- .summ(tt$P.Value, p.thresh, tt$logFC, logFC.thresh, dir)
		q <- .summ(tt$adj.P.Val, q.thresh, tt$logFC, logFC.thresh, dir)
		res[[dir]] <- rbind(p,q)
	}
	
	FC.thresh <- round(2^logFC.thresh,1)
	
	rn <- c()
	if(!is.null(p.thresh)) rn <- c(rn,paste0("P < ",   p.thresh))
	if(!is.null(q.thresh)) rn <- c(rn,paste0("FDR < ", q.thresh))
	cn <- paste0("FC > ", FC.thresh)
	res <- lapply(res, function(x) {dimnames(x) <- list(rn,cn); x})
	res
}


#' Run topTable on all contrasts from a linear model fit
#' Function to create a topTable for every contrast from an lmFit object, 
#' including the F-stat, if >1 contrasts
#' @param fit an \code{lmFit} object
#' @param number the number of DE genes. default=Inf which is all genes tested
#' @param genelist vector of genes to do the topTable on
#' @param adjust.method the multiple test correction method. default=\dQuote{BH}.
#'    see \code{\link[limma]{topTable}}
#' @param sort.by which column to sort by? see \code{\link[limma]{topTable}}
#' @param resort.by character string specifying statistic to sort the selected
#'        genes by in the output \code{data.frame}.  Possibilities are the
#'        same as for \code{sort.by}.
#' @param p.value a p.value threshold. default=1.0, ie all genes
#' @param lfc a logFC threhsold to use. default=0.0, ie all genes
#' @return a list of topTable objects, named by the coefficient names
#' @author Mark Cowley, 2010-10-07
#' @export
#' @importFrom limma topTable topTableF
topTable.all <- function(fit, number=Inf, genelist=fit$genes, adjust.method="BH",
         sort.by="p", resort.by=NULL, p.value=1.0, lfc=0.0) {
	res <- list()
	N <- ncol(fit$coefficients)
	for(i in 1:N) {
		res[[i]] <- topTable(fit, coef=i, number=number, genelist=genelist, adjust.method=adjust.method, sort.by=sort.by, resort.by=resort.by, p.value=p.value, lfc=lfc)
	}
	names(res) <- colnames(fit$coefficients)
	if( N > 1 ) {
		res$Fstat <- topTableF(fit, number=number, genelist=genelist, adjust.method=adjust.method, sort.by="F", p.value=p.value)
	}
	
	res
}



#' Is the argument a single topTable, or a list of topTables
#' @param tt a \code{data.frame}, or list of \code{data.frame}'s
#' @return logical: \code{TRUE} if \code{tt} is a list of \code{topTable}'s, \code{FALSE} otherwise
#' @author Mark Cowley, 2009-01-28
#' @export
is.topTable.list <- function(tt) {
	is.list(tt) && !is.data.frame(tt)
}




#' Add a moderated SE column to a topTable object.
#' 
#' moderated SE's are fit$stdev.unscaled * sqrt(fit$s2.post)
#' unmoderated SE's are fit$stdev.unscaled * fit$sigma
#' from p53 & 54 limma userguide:
#' for gene j, contrast k:
#' ujk = fit$stdev.unscaled = unscaled stdev
#' sj = fit$sigma = gene j's residual variance, with dj df
#' Bjk = linear model estimates = the logFC's
#' tstat.ord <- signal/SE <- Bjk/(ujk.sj) <- fit$coef/(fit$stdev.unscaled *
#' fit$sigma)
#' tstat.mod <- signal/SE' <- Bjk/(ujk.s'j) <- fit$coef/(fit$stdev.unscaled *
#' sqrt(fit$s2.post))
#' ^^^^ which follows a t-dist with d0 + dj d.f.
#' s2.post is the weighted avg of s2.prior and sigma^2 with weights
#' proportional to df.prior and df.residual, respectively.
#' 
#' @param fit see topTable
#' @param coef the coefficient index. if \code{NULL}, this is set to 1.
#' @param \dots arguments passed to \code{\link[limma]{topTable}}
#' @return a topTable object produced by \code{\link[limma]{topTable}}, in addition to an \dQuote{SE}
#'   column containing the moderated standard errors.
#' @author Mark Cowley, 2010-11-12
#' @export
#' @importFrom limma topTable
#' 
topTable.SE <- function(fit, coef=NULL, ...) {

	tt <- topTable(fit=fit, coef=coef, ...)
	if( is.null(coef) || (ncol(fit) == 1) )
		coef <- 1 

	if( !is.null(coef) ) {
		if( "s2.post" %in% names(fit) ) {
			SE <- fit$stdev.unscaled * sqrt(fit$s2.post)
		}
		else { # this will never be tested, since topTable only works on eBayes results.
			SE <- fit$stdev.unscaled * fit$sigma
		}
		if( is.matrix(SE) ) {
			SE <- SE[, coef]
		}
		SE <- SE[as.numeric(rownames(tt))]
		tt$SE <- SE
	}
	tt
}

#' Get the DE genes from a topTable.
#' Get the lists of DE genes that pass various P or Q thresholds, see also
#' \code{\link{summarise.topTable}}
#' 
#' @param tt a top.table, or a list of top tables
#' @param p.thresh vectors of thresholds for p.values
#' @param q.thresh vectors of thresholds for FDRs
#' @param lfc.thresh vector of POSITIVE log base 2 FC thresholds. Note that
#'   absolute logFC's are used
#' @param values also return the \sQuote{p}/\sQuote{q}/\sQuote{logFC} value?
#' @param as.df logical if \code{FALSE} the result will be a list, otherwise will be a
#'   \code{data.frame}
#' @param direction Which direction can the gene change be? one of \dQuote{either}, \dQuote{up} or \dQuote{down}
#' @return A list of DE genes from a toptable
#' @author Mark Cowley, 25/2/2008
#' @export
#' @seealso \code{\link{summarise.topTable}}
DEgenes.topTable <- function(tt, p.thresh=c(0.05, 0.001, 0.0001), 
						   q.thresh=c(0.25, 0.10, 0.05), 
						   lfc.thresh=c(0.585, 1, 2),
						   values=FALSE, as.df=FALSE,
						direction=c("either", "up", "down")) {
	
	direction <- direction[1]
	
	if( is.topTable.list(tt) ) {
		res <- lapply(tt, DEgenes.topTable,
					p.thresh=p.thresh, q.thresh=q.thresh, lfc.thresh=lfc.thresh, values=values, as.df=as.df, direction=direction)
		return( res )
	}
	# else if( length(direction) > 1 ) {
	# 	res <- list()
	# 	for(i in 1:length(direction)) {
	# 		res[[i]] <- DEgenes.topTable(tt=tt,p.thresh=p.thresh, q.thresh=q.thresh, lfc.thresh=lfc.thresh, values=values, as.df=as.df, direction=direction[i])
	# 	}
	# 	names(res) <- direction
	# 	return( res )
	# }
	else {
		if( direction != "either" && "logFC" %in% colnames(tt) ) {
			if( direction == "up" )
				tt <- tt[tt$logFC >= 0,]
			else if( direction == "down" )
				tt <- tt[tt$logFC < 0,]
			else
				stop("direction must be one of: either, up or down.\n")
		}
		
		IDcol <- which(colnames(tt) %in% c("ID", "ProbeSetID"))
		if( ! "logFC" %in% colnames(tt) )
			lfc.thresh <- numeric(0)

		ids <- list()
		if( length(p.thresh) > 0 ) {
			for(j in 1:length(p.thresh)) {
				ids[[j]] <- tt[ tt$P.Value < p.thresh[j], IDcol]
			}	
		}
	
		if( length(q.thresh) > 0 ) {
			qcol <- which( colnames(tt) %in% c("q", "q (Storey)") )
			if( length(qcol) == 0 )
				qcol <- which(colnames(tt) == "adj.P.Val")
			for(j in 1:length(q.thresh)) {
				ids[[j+length(p.thresh)]] <- tt[tt[, qcol] < q.thresh[j], IDcol]
			}
		}
	
		if( length(lfc.thresh) > 0 ) {
			for(j in 1:length(lfc.thresh)) {
				ids[[j+length(p.thresh)+length(q.thresh)]] <- tt[tt$logFC > lfc.thresh[j], IDcol]
			}
		}
	
		fc.thresh <- round(2^lfc.thresh,1)
		
		NAMES <- NA
		if( length(p.thresh) > 0 )
			NAMES <- c(NAMES, paste0("P<", p.thresh))
		if( length(q.thresh) > 0 )
			NAMES <- c(NAMES, paste0("q<", q.thresh))
		if( length(lfc.thresh) > 0 )
			NAMES <- c(NAMES, paste0("FC>", fc.thresh))
		NAMES <- na.rm(NAMES)
	
		names(ids) <- NAMES

		if( values ) {
			vals <- list()
			if( length(p.thresh) > 0 ) {
				for(j in 1:length(p.thresh)) {
					vals[[j]] <- tt$P.Value[ tt$P.Value < p.thresh[j] ]
				}
			}
			if( length(q.thresh) > 0 ) {
				for(j in 1:length(q.thresh)) {
					vals[[j+length(p.thresh)]] <- tt[ tt[,qcol] < q.thresh[j], qcol ]
				}
			}
			if( length(lfc.thresh) > 0 ) {
				for(j in 1:length(lfc.thresh)) {
					vals[[j+length(p.thresh)+length(q.thresh)]] <- tt$logFC[ tt$logFC > lfc.thresh[j] ]
				}
			}
			names(vals) <- NAMES
		
			res <- list()
			res[seq(1, length(ids)*2, 2)] <- ids
			res[seq(2, length(ids)*2, 2)] <- vals
			names(res) <- rep(NAMES, each=2)
		}
		else
			res <- ids

		if( as.df )
			res <- list2df(res)

		return( res )
	}
}


#' Export the DE genes from a topTable.
#' Get the DE genes passing various P or Q value threhsolds, and export them to
#' a xls file.
#' 
#' @param tt a top.table, or a list of top tables
#' @param file the path to the output tsv file
#' @param p.thresh vectors of thresholds for p.values
#' @param q.thresh vectors of thresholds for FDRs
#' @param lfc.thresh vector of POSITIVE log base 2 FC thresholds. Note that
#'   absolute logFC's are used
#' @param values also return the \sQuote{p}/\sQuote{q}/\sQuote{logFC} value?
#' @seealso \code{\link{DEgenes.topTable}}
#' 
#' @author Mark Cowley, 25/2/2008
#' @export
export.DEgenes.topTable <- function(tt, file=NULL, 
						p.thresh=c(0.05, 0.001, 0.0001), 
						q.thresh=c(0.25, 0.10, 0.05), 
						lfc.thresh=c(0.585, 1, 2),
						values=FALSE) {
	
	degenes <- DEgenes.topTable(tt, p.thresh=p.thresh, q.thresh=q.thresh, lfc.thresh=lfc.thresh, values=values, as.df=TRUE)
	
	write.xls(degenes, file, row.names=FALSE, na="")
}




#' Venn Diagram from 2 topTable objects
#' Plot a Venn Diagram of the differentially expressed genes within 2 topTable objects.
#' Generates a 4x3 panel of plots, 1 per each p.threshold, q.thresh, logFC.thresh and
#' top N sizes, where N defaults to 50, 100, 200.
#' 
#' @param tt1 a \code{data.frame} from \code{\link[limma]{topTable}} 
#' @param tt2 a \code{data.frame} from \code{\link[limma]{topTable}} 
#' @param p.thresh a vector of P Value thresholds, 1 per Venn Diagram
#' @param q.thresh a vector of FDR thresholds, 1 per Venn Diagram
#' @param logFC.thresh a vector of logFC thresholds, 1 per Venn Diagram
#' @param sizes a vector of top N sizes
#' @param names the names of the 2 topTable objects. default=\dQuote{A}, \dQuote{B}
#' @return none.
#' @export
#' @author Mark Cowley, 2011-08-02
#'
plot_venn2D_topTable <- function(tt1, tt2,  
	p.thresh=c(0.05, 0.001, 0.0001), 
	q.thresh=c(0.25, 0.10, 0.05), 
	logFC.thresh=c(0.585, 1, 2), 
	sizes=c(50, 100, 250), names=LETTERS[1:2]) {

	pop <- unionN(tt1$ID, tt2$ID)

	par(mfrow=c(4,3), cex.main=2)
	
	for(thresh in p.thresh)
		plot.venn2D(tt1$ID[tt1$P.Value < thresh], tt2$ID[tt2$P.Value < thresh], names=names, main=paste("P <", thresh), population=pop )

	for(thresh in q.thresh)
		plot.venn2D(tt1$ID[tt1$adj.P.Val < thresh], tt2$ID[tt2$adj.P.Val < thresh], names=names, main=paste("q <", thresh), population=pop )
	
	for(thresh in logFC.thresh) {
		FC.thresh <- round(2^thresh, 1)
		plot.venn2D(tt1$ID[abs(tt1$logFC) > thresh], tt2$ID[abs(tt2$logFC) > thresh], names=names, main=paste("FC >", FC.thresh), population=pop )
	}

	for(N in sizes)
		plot.venn2D(tt1$ID[1:N], tt2$ID[1:N], names=names, main=paste("top", N), population=pop )
}



#' Venn Diagram from 3 topTable objects
#' Plot a Venn Diagram of the differentially expressed genes within 3 topTable objects.
#' Generates a 4x3 panel of plots, 1 per each p.threshold, q.thresh, logFC.thresh and
#' top N sizes, where N defaults to 50, 100, 200.
#' 
#' @param tt1 a \code{data.frame} from \code{\link[limma]{topTable}}
#' @param tt2 a \code{data.frame} from \code{\link[limma]{topTable}}
#' @param tt3 a \code{data.frame} from \code{\link[limma]{topTable}}
#' @param p.thresh a vector of P Value thresholds, 1 per Venn Diagram
#' @param q.thresh a vector of FDR thresholds, 1 per Venn Diagram
#' @param logFC.thresh a vector of logFC thresholds, 1 per Venn Diagram
#' @param sizes a vector of top N sizes
#' @param names the names of the 3 topTable objects. default=\dQuote{A}, \dQuote{B}, \dQuote{C}
#' @return none.
#' @export
#' @author Mark Cowley, 2011-08-02
plot_venn3D_topTable <- function(tt1, tt2, tt3, p.thresh=c(0.05, 0.001, 0.0001), q.thresh=c(0.25, 0.10, 0.05), logFC.thresh=c(0.585, 1, 2), sizes=c(50, 100, 250), names=LETTERS[1:3]) {

	pop <- unionN(tt1$ID, tt2$ID, tt3$ID)
	
	par(mfrow=c(4,3), cex.main=2)
	
	for(thresh in p.thresh)
		plot.venn3D(tt1$ID[tt1$P.Value < thresh], tt2$ID[tt2$P.Value < thresh], tt3$ID[tt3$P.Value < thresh], names=names, main=paste("P <", thresh), population=pop )
		
	for(thresh in q.thresh)
		plot.venn3D(tt1$ID[tt1$adj.P.Val < thresh], tt2$ID[tt2$adj.P.Val < thresh], tt3$ID[tt3$adj.P.Val < thresh], names=names, main=paste("q <", thresh), population=pop )
	
	for(thresh in logFC.thresh) {
		FC.thresh <- round(2^thresh, 1)
		plot.venn3D(tt1$ID[abs(tt1$logFC) > thresh], 
					tt2$ID[abs(tt2$logFC) > thresh], 
					tt3$ID[abs(tt3$logFC) > thresh], names=names, main=paste("FC >", FC.thresh), population=pop )
	}

	for(N in sizes)
		plot.venn3D(tt1$ID[1:N], tt2$ID[1:N], tt3$ID[1:N], names=names, main=paste("top", N), population=pop )
}


#' A topTable volcano plot.
#' Volcano plot of a toptable - similar to \code{\link[limma]{volcanoplot}} from \code{limma},
#'  but more flexible in terms of what is plotted on the y-axis.
#' 
#' @param tt a \code{data.frame} object from calling \code{\link[limma]{topTable}} with all genes.
#' @param yaxis choose one of \dQuote{p}, \dQuote{q}, \dQuote{B}, \dQuote{t}, 
#' \dQuote{absT} to plot the raw P-values, the q
#'   value (either from a column called \dQuote{q}, or \dQuote{adj.P.Val} 
#' in that order), the B statistics (log odds), the t-statistic, or the absolute 
#' t-statistic, respectively.
#' @param xaxis the x-axis type. one of: \dQuote{logFC} (default), \dQuote{signedFC}
#'    (see \code{\link{logFC2signedFC}}), or \dQuote{FC} (unlogged Fold Change)
#' @param yThresh For colouring significant points. For eg, choose 0.05 for
#'   P/q, or 3 for B, or 5 for t
#' @param lfc the absolute log fold change threshold. eg 0.585 for 1.5 FC, or
#'   1.0 for a 2-fold change
#' @param highlight the number of points to highlight by name
#' @param names the names to use if highlight > 0. This should have as many
#'   rows as were passed into \code{lmFit}, and should be in the same order as the
#'   genes passed into lmFit. The tt's row names (which are numeric) will be
#'   used to index into this vector of names.
#' @param main the plot title.
#' @param cex.points see \code{\link{par}}
#' @param ablines whether to add vertical dashed lines at \code{yThresh} and \code{xThresh}.
#'    \dQuote{}, \dQuote{x}, \dQuote{y}, or \dQuote{xy} for none, x (ie vertical lines at +/- xThresh only), 
#'    y (ie horizontal line at yThresh only), x and y
#' @param abline.col the colour of the ablines, if \code{ablines != ""}
#' @param colour.scheme Control the higlighting of DE genes, via either black/red dots (\dQuote{red}) 
#'    or closed/open circles(\dQuote{black})
#' @param xlim see \code{\link{par}} if \code{NULL} it will be made symmetrical so that all data points fit
#' @param ylim see \code{\link{par}}. if \code{NULL} it will fit all data points
#' @param xlab The x-axis label. if \code{NULL} it's auto-determined from \code{xaxis}
#' @param ylab The y-axis label. if \code{NULL} it's auto-determined from \code{yaxis}
#' @param do.par logical: setup the plot parameters
#' @param \dots further arguments passed to \code{\link{plot}}
#' @return nothing
#' @author Mark Cowley, 2008-07-25
#' @export
volcanoplot_topTable <- function(tt, 
		yaxis=c("p","q","B", "t"),
		xaxis=c("logFC", "signedFC", "FC"),
		yThresh=0.05, lfc=log2(1.5), 
		highlight=0, names=tt$ID, 
		main="Volcano Plot", 
		cex.points=0.2, 
		ablines="", 
		abline.col="grey",
		colour.scheme="red",
		xlim=NULL, ylim=NULL, 
		xlab=NULL, ylab=NULL,
		do.par=TRUE, ...) {
	
	if( ! "logFC" %in% colnames(tt) )
		stop("logFC column not in this tt; are you trying to volcano plot a topTable from an Ftest?")

	y <- rep(0, nrow(tt))
	YLAB <- ylab # IF ylab was set as an argument, then remember it's value and over-write it after this code block (far fewer if is.null statements this way).
	ylab <- "ERROR"
	yaxis <- yaxis[1]
	if( yaxis %in% c("p", "P", "P.Value") ) {
		y <- -log10(tt$P.Value)
		ylab="P-value (-log10)"
		yThresh <- -log10(yThresh)
	}
	else if( yaxis %in% c("q", "Q", "FDR") ) {
		if( "q" %in% colnames(tt) )
			y <- tt$q
		else
			y <- tt$adj.P.Val
		y <- -log10(y)
		ylab="FDR"
		yThresh <- -log10(yThresh)
	}
	else if( yaxis %in% c("B") ) {
		y <- tt$B
		ylab="log odds"
	}
	else if( yaxis %in% c("t", "TRUE") ) {
		y <- tt$t
		ylab="t statistic"
	}
	else if( yaxis %in% c("absT", "abst") ) {
		y <- abs(tt$t)
		ylab="|t statistic|"
	}
	else {
		stop("invalid yaxis; try one of P, q or B")
	}
	if( !is.null(YLAB) )
		ylab <- YLAB
	
	x <- rep(0, nrow(tt))
	XLAB <- xlab # IF ylab was set as an argument, then remember it's value and over-write it after this code block (far fewer if is.null statements this way).
	
	xlab <- "ERROR"
	xaxis <- xaxis[1]
	if( xaxis %in% "logFC" ) {
		xlab <- "log2 Fold Change"
		x <- tt$logFC
	}
	else if( xaxis %in% "signedFC" ) {
		xlab <- "Signed Fold Change"
		x <- logFC2signedFC(tt$logFC)
	}
	else if( xaxis %in% "FC" ) {
		xlab <- "Fold Change"
		x <- 2^tt$logFC
	}
	else if( xaxis %in% "absFC" ) {
		xlab <- "Fold Change"
		x <- 2^abs(tt$logFC)
	}
	else {
		stop("invalid xaxis; try one of logFC, signedFC, FC")		
	}
	if( !is.null(XLAB) )
		xlab <- XLAB
	
	non.na.idx <- !is.na(x) & !is.na(y)
	x <- x[non.na.idx]
	y <- y[non.na.idx]
	tt <- tt[non.na.idx, ]
	
	# colour points that are more extreme than the given thresholds
	if( "col" %in% colnames(tt) ) {
		col <- tt$col
	}
	else {
		col <- rep(1, nrow(tt))
		pch <- rep(16, nrow(tt))
		
		yPaint <- NULL
		xPaint <- NULL
		if( !is.null(yThresh) ) {
			# if( ! yaxis %in% c("t", "TRUE") )
				yPaint <- which(abs(y)>yThresh)
			# else
			# 	yPaint <- which(abs(y)>yThresh)
			
		}
		if( !is.null(lfc) )
			xPaint <- which(abs(x)>lfc)
		if( colour.scheme == "black" ) {
			pch <- rep(16, nrow(tt))
			pch[intersect(yPaint, xPaint)] <- 1
		}
		else {
			col[intersect(yPaint, xPaint)] <- "red"
		}
		# col[intersect(yPaint, xPaint)] <- 2
	}
	
	opar <- NULL
	if( do.par ) {
		opar <- par(c("las", "cex.axis", "cex.lab", "cex.main", "mar"))
		par(las=1, cex.axis=1.5, cex.lab=1.5, cex.main=1.6, mar=par()$mar+1)
	}

	if( is.null(ylim) ) {
		ylim <- range(y, na.rm=TRUE)
		if( yaxis %in% c("p", "P", "P.Value", "q", "Q", "FDR") ) {
			ylim <- c(0, max(ylim[2], 2))
		}
	}
		

	# # add a small red tick on the axis at the points where the thresholds have been set.
	# axis(side=1, at=c(-lfc, lfc), col="red", tick=TRUE, labels=FALSE)
	# axis(side=2, at=c(yThresh), col="red", tick=TRUE, labels=FALSE)
	
	if( is.null(xlim) ) {
		xlim <- symmetricise(range(x, na.rm=TRUE))
	}
	plot(x, y, main=main, 
		 xlab=xlab, ylab=ylab, 
		 col=col, xlim=xlim, ylim=ylim,
		 pch = pch, cex = cex.points, 
		 yaxt=ifelse(yaxis %in% c("t", "T"), "s","n"), ...)

	if( yaxis %in% c("t", "T") ) {
		abline(h=0, v=0, lty="dashed", col=abline.col)
	}
	else {
		yTicks <- sort(unique(c(yThresh, -log10(0.05), seq(0,ceiling(max(y)), 1))))
		yTicks <- setdiff(yTicks, 1)
		# cat(yTicks, "\n")
		axis(side=2, at=yTicks, labels=prettyNum(10^-yTicks), line=0)
	}
	
	if (highlight > 0) {
        if (is.null(names)) 
            names <- 1:length(x)
        names <- as.character(names)
        o <- order(y, decreasing = TRUE)
        i <- o[1:highlight]
        text(x[i], y[i], labels = substring(names[i], 1, 8), 
            cex = 0.8, col = "blue")
    }

	if( grepl("x", ablines) )
		abline(v=lfc*c(-1,1), lty="dashed", col=abline.col)
	if( grepl("y", ablines) )
		abline(h=yThresh, lty="dashed", col=abline.col)

	if( do.par ) {
		par(opar)
	}
}


#' A wrapper function to plot 2 volcano plots side by side, for the logFC vs P
#' and q values.
#' 
#' @param tt a \code{data.frame} object from calling \code{\link[limma]{topTable}} with all genes.
#' @param pThresh For colouring significant points in the p-value panel., default=0.001
#' @param qThresh For colouring significant points in the FDR panel., default=0.05
#' @param lfc the absolute log fold change threshold. eg 0.585 for 1.5 FC, or
#'   1.0 for a 2-fold change
#' @param cex.points see \code{\link{par}}
#' @param ablines whether to add vertical dashed lines at \code{yThresh} and \code{xThresh}.
#'    \dQuote{}, \dQuote{x}, \dQuote{y}, or \dQuote{xy} for none, x (ie vertical lines at +/- xThresh only), 
#'    y (ie horizontal line at yThresh only), x and y
#' @param \dots further arguments passed to \code{\link{volcanoplot_topTable}}
#' @seealso \code{\link{volcanoplot_topTable}}
#' @author Mark Cowley, 2008-08-01
#' @export
volcanoplot_topTable_PandQ <- function(tt, pThresh=0.001, qThresh=0.05, lfc=log2(1.5), cex.points=0.5, ablines="", ...) {
	par(mfrow=c(1,2))
	volcanoplot_topTable(tt, "p", yThresh=pThresh, lfc=lfc, cex.points=cex.points, ablines=ablines, ...)
	volcanoplot_topTable(tt, "q", yThresh=qThresh, lfc=lfc, cex.points=cex.points, ablines=ablines, ...)
}


#' Convert logFC (FoldChange) to signed FC
#' When logFC<0, unlogged FC values get squashed in (0,1). The signed FC
#' has the same magnitude whether the logFC was >0 or <0, but if <0, a leading - is
#' added
#' @param fc a numeric vector of log-base-2 Fold Change values
#' @return a numeric vector of signed fold change values
#' @author Mark Cowley, 2009-02-04
#' @export
logFC2signedFC <- function(fc) {
	res <- rep(NA, length(fc))
	res[fc >= 0] <- 2^(fc[fc >= 0])
	res[fc < 0] <- -1 * (2^abs(fc[fc < 0]))
	res[fc == 0] <- 1
	res
}




#' 'fix' the logFC column for biologists (sorry!)
#' 
#' @param tt the output from \code{\link[limma]{topTable}}
#' @param digits how many digits? see \code{\link{round}}
#' @return an edited topTable, with 3 extra columns (eg for logFC=-1.82):\cr
#' \dQuote{FC}: the unlogged fold change (eg 0.283)\cr
#' \dQuote{direction}: up or down (eg down)\cr
#' \dQuote{absFC}: the unlogged absolute fold change (eg 3.531)\cr
#'
#' @author Mark Cowley, 2008-07-14
#' @examples
#' \dontrun{
#' head( topTable.fixFC(tt) )
#' #            ID logFC     t  P.Value adj.P.Val     B    FC direction absFC
#' # 16143 7923547  2.58 13.44 3.31e-06    0.0695 -1.41 5.985        up  5.99
#' # 12361 8096301  2.59  9.52 3.22e-05    0.3382 -1.62 6.039        up  6.04
#' # 16282 7953200  1.74  8.62 6.11e-05    0.4273 -1.70 3.331        up  3.33
#' # 540   7995729  1.64  8.03 9.52e-05    0.4996 -1.77 3.109        up  3.11
#' # 8641  8141374 -1.42 -7.43 1.55e-04    0.6164 -1.85 0.374      down  2.67
#' # 11546 8007848  1.48  7.28 1.76e-04    0.6164 -1.87 2.790        up  2.79
#' }
#' @export
topTable.fixFC <- function(tt, digits=4) {
	if( "logFC" %in% colnames(tt) ) {
		# tt$FC <- round(2^tt$logFC, digits)
		tt$absFC <- round(2^abs(tt$logFC), digits)
		tt$direction <- ifelse(tt$logFC >= 0, "up", "down")
		tt <- move.column(tt, c("absFC", "direction"), which(colnames(tt) == "logFC")+1)
	}
	tt
}




#' Export a topTable, or a list of top tables to an excel file.
#' Export a topTable, or a list of top tables to an excel file, with an
#' optional summary of the number of DE genes passing various generic
#' thresholds.
#' 
#' The first worksheet will be the summary, and subsequent worksheets are for
#' each top table supplied.
#' Additional benefits:\cr
#' - Limits the number of significant figures to digits correctly for the
#' P and q values\cr
#' - includes a Pcount column\cr
#' - includes gene information which uses the ID column and rownames of annot
#' object\cr
#' - includes a few columns of coefficients, for eg if the group-means
#' parameterisation
#' is used, then the values of the coefficients (ie the value for each 'group')
#' may be appropriate to export.
#' 
#' @param tt a topTable data.frame
#' @param file the output file name
#' @param annot the probe-level annotation object. rownames should be set to
#'   the probe ID, so tt$ID can match them
#' @param fixFC fix the logFC to produce an absFC and direction columns
#' @param Pcount an optional vector of Pcounts (ie how many times each probe
#'   was detected as present). Should be named so that tt$ID can match the
#'   probes
#' @param coefficients an optional data.frame, of coefficients from the lmFit.
#'   You can use this arg to add custom columns to the topTable. eg the AvgExpr
#' @param fit the lmFit object. Used to obtain the Standard error of each
#'   probe. using fit$stdev.unscaled * fit$sigma
#' @param digits how many digits to round to?
#' @param summary insert a worksheet into the excel workbook which contains a
#'   summary of the toptable
#' @param drop.Bstat logical. If TRUE then the "B" column is not exported.
#' @param adj.Pval.colname rather than "adj.P.Val", you can rename to column to
#'   "FDR", or "q", or similar.
#' @param \dots additional arguments passed to write.xls
#' @author Mark Cowley
#' @export
#' @importFrom excelIO write.xls
export.topTable <- function(tt, file, annot=NULL, fixFC=TRUE, Pcount=NULL, coefficients=NULL, fit=NULL, digits=4, summary=TRUE, drop.Bstat=TRUE, adj.Pval.colname="FDR", ... ) {
	# first calculate a summary of the toptable
	if( summary ) {
		smry <- summarise.topTable(tt)
	}

	if( is.topTable.list(tt) ) {
		#
		# then re-run export.topTable on all elements of the list, but do not write
		# them out to individual files -- thus file=NULL
		# The list of parsed/edited topTables will be exported once after the following
		# large else{} block.
		#
		tt <- lapply(tt, export.topTable, file=NULL,
			annot=annot, fixFC=fixFC, Pcount=Pcount, coefficients=coefficients, fit=fit, digits=digits, summary=FALSE, drop.Bstat=drop.Bstat, adj.Pval.colname=adj.Pval.colname, ... )
		if( summary ) {
			smry <- list(smry)
			names(smry) <- "summary"
			tt <- lbind(smry, tt)
		}
	}
	else {
		if( fixFC )
			tt <- topTable.fixFC(tt, digits)

		tt <- round.topTable(tt, digits)

		if( !is.null(Pcount) ) {
			if( is.vector(Pcount) ) {
				Pcount <- as.data.frame(Pcount)
			}
			tt <- merge(tt, Pcount, by.x="ID", by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)

			idx <- min(which(colnames(tt) %in% c("t", "F")))
			tt <- move.column(tt, colnames(Pcount), idx)
		}
		
		#
		# insert the coefficients from the model fit -- rather than AveExpr which
		# is a bit dull...
		#
		if( !is.null(coefficients) ) {
			if( length(dim(coefficients)) == 0 ) { # ie, coef is a vector
				if( length(coefficients) == nrow(tt) ) {
					coefficients <- as.matrix(coefficients)
					colnames(coefficients) <- "coef"
				}
				else {
					stop("Supplied coefficients are not a matrix, and don't have the same length as the toptable.\n")
				}
			}

			if( is.null(colnames(coefficients)) )
				colnames(coefficients) <- "coeff"
			coefficients <- round(coefficients, digits)

			tt <- merge(tt, coefficients, 
						by.x="ID", by.y="row.names",
						all.x=TRUE, all.y=FALSE, sort=FALSE)

			idx <- min(which(colnames(tt) %in% c("AveExpr", "t", "F")))
			tt <- move.column(tt, colnames(coefficients), idx)
		}

		if( !is.null(fit) ) {
			# then add an SE column.
			err <- fit$stdev.unscaled * fit$sigma
			err <- round(err, digits)
			if( is.vector(err) ) {
				err <- data.frame(SE=err)
				rownames(err) <- rownames(fit$stdev.unscaled)
			}
			else if( ncol(err) == 1 ) {
				colnames(err) <- "SE"
			}
			else {
				colnames(err) <- paste(colnames(err), "SE", sep=".")
			}
			
			tt <- merge(tt, err, 
						by.x="ID", by.y="row.names",
						all.x=TRUE, all.y=FALSE, sort=FALSE)
			# tt <- cbind(tt, err[match(tt$ID, rownames(err)), ])
			idx <- min(which(colnames(tt) %in% c("AveExpr", "t", "F")))
			tt <- move.column(tt, colnames(err), idx)
		}

		if( !is.null(annot) ) {
			tt <- merge(tt, annot, 
					by.x="ID", by.y=colnames(annot)[1], 
					all.x=TRUE, all.y=FALSE, sort=FALSE)
		}

		if( colnames(tt)[1] == "ID" )
			colnames(tt)[1] <- "ProbeSetID"

		if( !is.null(adj.Pval.colname) && "adj.P.Val" %in% colnames(tt) )
			colnames(tt)[colnames(tt) == "adj.P.Val"] <- adj.Pval.colname

		if( drop.Bstat && "B" %in% colnames(tt) ) {
			tt$B <- NULL
		}

		if( summary )
			tt <- list(summary=smry, topTable=tt)
	}

	#
	# if processing a toptable list, then we don't want to write a file until the end
	# thus file has been set to NULL, and the invisible result will be cauught by the 
	# lapply at the start of this function.
	#
	if( !is.null(file) ) {
		write.xls(tt, file, row.names=FALSE, ...)
	}
	
	invisible(tt)
}
# CHANGELOG
# 11/3/2008: v1
# 2009-01-28: extensively updated to export xls files, and handle lists of toptables
# 2010-01-19: support added for 'fit' which inserts SE's into the output.


#' Round the numeric columns in a topTable
#' Take a top table, and round the columns to \code{digits} places
#' 
#' @param x a topTable object
#' @param digits the number of digits to round to
#' @return a topTable object with numerical columns rounded.
#' @author Mark Cowley
#' @method round topTable
round.topTable <- function(x, digits=4) {
	if( "F" %in% colnames(x) ) {
		# identify the numeric columns to the left of, and upto "F" which contain the contrasts
		numeric.cols <- colclasses(x) == "numeric"
		Fcol <- which(colnames(x)=="F")
		left.of.F <- (1:ncol(x) > 1) & (1:ncol(x) <= Fcol)
		cols <- which( numeric.cols & left.of.F )
		# cols <- intersect(2:which(colnames(x)=="F"), which(colapply(x, is.numeric)))
		for(col in cols)
			x[,col] <- round(x[,col], digits)
	}
	else if( "t" %in% colnames(x) ) {
		x$logFC <- round(x$logFC, digits)
		x$t <- round(x$t, digits)
	}
	if("AveExpr" %in% colnames(x)) x$AveExpr <- round(x$AveExpr, digits)
	if("B" %in% colnames(x)) x$B <- round(x$B, digits)
	if("SE" %in% colnames(x)) x$SE <- round(x$SE, digits)

	if("P.Value" %in% colnames(x)) x$P.Value <- format.pval(x$P.Value, digits)
	if("adj.P.Val" %in% colnames(x)) x$adj.P.Val <- format.pval(x$adj.P.Val, digits)
	if("q" %in% colnames(x)) x$q <- format.pval(x$q, digits)
	if("FDR" %in% colnames(x)) x$FDR <- format.pval(x$FDR, digits)

	return(x)
}
# CHANGELOG
# 2012-07-06: added @S3method
# 2012-07-24: @method is what works

#' Take a toptable (Fstat or t-stat) and make a label for each row.
#' 
#' @param tt a toptable
#' @param probes if \code{NULL}, all probes from \code{tt$ID} are used. otherwise, probes can
#'   be a numeric vector (eg 1 10) or character vector corresponding to
#'   \code{tt$ID}
#' @param probe2genesymbol a 2 column \code{data.frame} with probe ID's in column 1, and
#'    gene symbols in column 2. If there's no gene symbol for a probe, then there should 
#'    be an \code{NA}.
#' @return a \code{character vector} with same length as probes, eg: 
#' [1] "probeset:10543233 F=1964 P=7.17e-13 FDR=1.16e-10"\cr
#' [2] "probeset: 10410211 t=9.26 P=8.49e-06 FDR=0.0193"\cr
#' @author Mark Cowley, 2009-07-16
#' @export
topTable.to.label <- function(tt, probes=NULL, probe2genesymbol=NULL) {
	if( is.null(probes) ) {
		probes <- tt$ID
	}
	else if( all(is.character(probes)) ) {
		tt <- tt[match(probes, tt$ID), ]
	}
	else if( all(is.numeric(probes)) )
		tt <- tt[probes,]
	
	if( !is.null(probe2genesymbol) ) {
		probe2genesymbol <- probe2genesymbol[match(probes, probe2genesymbol[,1]),2]
		probe2genesymbol[is.na(probe2genesymbol)] <- "---"
		if("F" %in% colnames(tt)) {
			paste(probe2genesymbol, " (", tt$ID, ")\n",
					" F=", prettyNum(tt$F),
					" P=", prettyNum(tt$P.Value),
					" FDR=", prettyNum(tt$adj.P.Val), sep="")
		}
		else if("t" %in% colnames(tt)) {
			paste(probe2genesymbol, " (", tt$ID, ")\n",
					" t=", prettyNum(tt$t),
					" P=", prettyNum(tt$P.Value),
					" FDR=", prettyNum(tt$adj.P.Val), sep="")
		}
	}
	else {
		if("F" %in% colnames(tt)) {
			paste("probeset: ", tt$ID, 
					"\n",
					" F=", prettyNum(tt$F),
					" P=", prettyNum(tt$P.Value),
					" FDR=", prettyNum(tt$adj.P.Val), sep="")
		}
		else if("t" %in% colnames(tt)) {
			paste("probeset: ", tt$ID, 
					"\n",
					" t=", prettyNum(tt$t),
					" P=", prettyNum(tt$P.Value),
					" FDR=", prettyNum(tt$adj.P.Val), sep="")
		}
	}
}

