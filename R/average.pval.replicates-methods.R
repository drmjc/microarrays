#' average.pval.replicates
#' 
#' @details Given a 2D table of pvalues, such as the detection pvalue for each probe for
#' each sample, compute an average p-value for those samples which are replicates
#' of each other.\cr
#' This uses the more correct calculation of
#' \code{p*(1-ln(p)), where p = p1*p2}, which is easily generalised to >2 replicates.
#'
#' @note Since detection pvals=0.00000 are quite common, and log(0) == NaN,
#' if \emph{any} of the p-values are 0, then the averaged p-value is also set to 0. I'm not
#' 100\% sure this is the best approach -- it's certainly correct of all values are 0, but
#' theoretically, 1 sample could have P=1.0 and the other P=0.
#' 
#' @param x a \code{matrix} or \code{data.frame} of p-values in [0,1]
#' @param classes a \code{character}, \code{numeric} or \code{factor} vector of classes that 
#' each sample belongs to. in this instance, class usually represents the original
#' sample from which each array has been generated from.
#' @return a \code{matrix} or \code{data.frame} of 'averaged' p.values. See Details.
#' 
#' @author Mark Cowley, Mark Pinese
#' 
#' @exportMethod average.pval.replicates
#' @rdname average.pval.replicates-methods
#' 
#' @docType methods
setGeneric(
	"average.pval.replicates",
	function(x, classes) {
		standardGeneric("average.pval.replicates")
	}
)

#' @rdname average.pval.replicates-methods
#' @aliases average.pval.replicates,matrix,factor-method
setMethod(
	"average.pval.replicates",
	signature=signature("matrix", "factor"),
	function(x, classes) {
		length(classes) == ncol(x) || stop("length(classes) != ncol(x)")
		res <- matrix(NA, nrow(x), length(levels(classes)), dimnames=list(rownames(x), levels(classes)))
		# there must be a convenient way to use by/apply/ddply or similar to do this on wide matrices...
		for(i in seq(along=levels(classes))) {
			class <- levels(classes)[i]
			if( sum(classes == class) == 1 ) {
				p <- x[, which(classes==class)]
			}
			else {
				p <- apply(x[, classes==class], 1, prod)
				p <- p*(1-log(p))
				p[is.nan(p)] <- 0 # since log(0)==NaN
			}
			res[, i] <- p
		}
		
		res
	}
)
# eg output from averaging the 2 replicates:
#              APGI_2222-TR.1 APGI_2222-TR.2      avg
# ILMN_1762337         0.0000         0.0000 0.000000
# ILMN_2055271         0.0506         0.0182 0.007358
# ILMN_1736007         0.3558         0.4130 0.428768
# ILMN_2383229         0.0000         0.0000 0.000000
# ILMN_1779670         0.5195         0.6922 0.727374
# ILMN_1653355         0.1260         0.1130 0.074756
# ILMN_1717783         0.7701         0.8208 0.922052
# ILMN_1705025         0.2156         0.1312 0.129108
# ILMN_1814316         0.2494         0.1234 0.137866
# ILMN_2359168         0.7494         0.7935 0.903720
# ILMN_1731507         0.8649         0.9247 0.978467
# ILMN_1787689         0.8390         0.6455 0.873670
# ILMN_3241953         0.0000         0.0000 0.000000
# ILMN_2136495         0.6442         0.8662 0.883529
# ILMN_1668111         0.1481         0.4571 0.249942
# ILMN_2295559         0.2857         0.6896 0.517088
# ILMN_1735045         0.0000         0.0000 0.000000
# ILMN_1680754         0.0195         0.2935 0.035245
# ILMN_2375184         0.7156         0.5338 0.749568
# ILMN_1659452         0.0779         0.0416 0.021803
# ILMN_1767388         0.5987         0.1117 0.247748
# ILMN_1675204         0.8753         0.9143 0.978581
# ILMN_1673870         0.8429         0.7532 0.923318
# ILMN_1755321         0.0013         0.0000 0.000000
# ILMN_1698554         0.0039         0.0000 0.000000
# ILMN_1814092         0.0013         0.0000 0.000000
# ILMN_1760414         0.0013         0.0143 0.000221
# ILMN_2061446         0.0000         0.0000 0.000000
# ILMN_1752884         0.3117         0.3584 0.356588
# ILMN_1660703         0.6701         0.4481 0.661496
# ILMN_1668851         0.7455         0.8195 0.911958
# ILMN_2270015         0.5857         0.5351 0.677025
# ILMN_1809959         0.9974         0.8273 0.983731
# ILMN_2357031         0.7922         0.6844 0.874091


#' @rdname average.pval.replicates-methods
#' @aliases average.pval.replicates,matrix,character-method
setMethod(
	"average.pval.replicates",
	signature=signature("matrix", "character"),
	function(x, classes) {
		classes <- factor(classes, levels=unique(classes))
		average.pval.replicates(x,classes)
	}
)

#' @rdname average.pval.replicates-methods
#' @aliases average.pval.replicates,data.frame,ANY-method
setMethod(
	"average.pval.replicates",
	signature=signature("data.frame", "ANY"),
	function(x, classes) {
		res <- average.pval.replicates(as.matrix(x), classes)
		res <- as.data.frame(res)
		res
	}
)

#' @rdname average.pval.replicates-methods
#' @aliases average.pval.replicates,ANY,missing-method
setMethod(
	"average.pval.replicates",
	signature=signature("ANY", "missing"),
	function(x, classes) {
		stop("You must supply a classes vector, same length as ncol(x)")
	}
)
