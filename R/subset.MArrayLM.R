#' Subset an MArrayLM object.
#' S3 method to subset the rows of an MArrayLM object, a la \code{\link[base]{subset}}
#' 
#' @param x an MArrayLM object to be subsetted.
#' @param subset a logical vector of rows/observations to include. Different to
#'   the default subset, you can also specify a character vector, or a numeric
#'   vector of indices. Leave missing if you want to ignore this argument.
#' @param select a logical vector of columns/conditions to include. Like
#'   subset, you can specify a character or numeric argument. Don't try to set
#'   this if ncol(x) == 1l you'll get a warning and this arg will be ignored.
#' @param drop ignored.
#' @param \dots ignored.
#' @return an MArrayLM with less rows, and or columns.  
#' @seealso \code{\link{subset}} \code{\link{subset.data.frame}}
#' @author Mark Cowley, 2009-11-10
#' 
#' @importClassesFrom limma MArrayLM
#' @S3method subset MArrayLM
subset.MArrayLM <- function(x, subset, select, drop=FALSE, ...) {
	
	# convert subset into a logical
	if( missing(subset) ) {
		r <- TRUE
	}
	else if( is.character(subset) ) {
		if( is.data.frame(x$genes) )
			r <- x$genes$ID %in% subset
		else
			r <- x$genes %in% subset
	}
	else if( is.numeric(subset) && max(subset) < length(x$coef) && min(subset) > 0 ) {
		r <- rep(FALSE, length(x$sigma))
		r[subset] <- TRUE
	}
	else if( is.logical(subset) ) {
		r <- subset
	}
	else {
		stop("subset must either be an integer vector of indicies, or a character vector, or a logical vector.\n")
	}
	if( all(!r) )
		stop("You've excluded all measurements.\n")

	# convert select into a logical.
	if( missing(select) ) {
		vars <- TRUE
	}
	else if( ncol(x) == 1 ) { # thus, select HAS been specified.
		vars <- TRUE
		warning("Trying to select conditions when there is only one to begin with. Ignoring this argument.\n")
	}
	else if( is.character(select) ){
		vars <- rep(FALSE, ncol(x))
		vars[match(select, colnames(x))] <- TRUE
	}
	else if( is.numeric(select) && max(select) < ncol(x) && min(select) > 0 ) {
		vars <- rep(FALSE, ncol(x))
		vars[select] <- TRUE		
	}
	else if( is.logical(select) ) {
		vars <- select
	}
	else {
		stop("select must either be an integer vector of indicies, or a character vector, or a logical vector.\n")
	}
	if( all(!vars) )
		stop("You've excluded all conditions.\n")
	
	
	#
	# Do the subsetting.
	#
	if( ncol(x) > 1 ) {
		res <- x[r, vars] # this works just fine; Smyth et al must have written a "[<-" method
	}
	else {
		#
		# the $genes is quite tricky to get right. if you print an MArrayLM, it looks like $genes is just a vector.
		# It's actually a data.frame with names = "ID"...
		# However when I manipulate $genes above, when you print the resulting MArrayLM, then it clearly looks like
		# a data.frame with 1 column called ID...
		#
		if( ncol(x$genes) > 1 )
			genes <- x$genes[r,]
		else {
			genes <- x$genes$ID[r]
			genes <- as.data.frame(list(ID=genes))
		}
		
		res <- x
		res$coefficients <- res$coefficients[r]
		res$stdev.unscaled <- res$stdev.unscaled[r]
		res$sigma <- res$sigma[r]
		res$df.residual <- res$df.residual[r]
		res$genes <- genes
		res$s2.post <- res$s2.post[r]
		res$t <- res$t[r]
		res$p.value <- res$p.value[r]
		res$lods <- res$lods[r]
		res$F <- res$F[r]
		res$F.p.value <- res$F.p.value[r]
	}
	
	return( res )
}
