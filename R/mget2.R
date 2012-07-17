#' enhanced mget
#' 
#' Enchanced mget to look up x in y. Changes in addition to \code{\link{mget}}
#' include:\cr
#' * Values can be sorted alphanumerically\cr
#' * Since some keys have multiple values, this function returns all, or just
#' the first (hint; use \code{first=TRUE, sort=TRUE}), assuming that the first
#' is the oldest. \cr
#' * keys in \code{x} are allowed to be \code{NA} (hint: \code{na.rm=FALSE})
#' 
#' @section TODO: provide more sorting options to choose the newest.
#' 
#' @param x a vector of keys, \code{NA}'s allowed.
#' @param y an env, or AnnDbBimap, or ProbeAnnDbBimap
#' @param sort logical: if \code{TRUE}, then sort the values for each key
#' @param first logical: if \code{TRUE}, then choose the first value for 
#'  each key, after the optional sort
#' @param na.rm logical: remove those keys that had no matches? if \code{TRUE}, 
#'   a vector with no \code{NA}'s will be returned. if \code{FALSE}, then return
#' a list of key to value maps, possibly with \code{NA}'s.
#' @return either a vector or list of values, depending on the value of 
#'  \code{na.rm}.
#' @author Mark Cowley, 2011-08-25
#' 
#' @export
#' @importFrom AnnotationDbi mget
#' 
#' @examples
#' if( require(org.Hs.eg.db) ) {
#'   sym <- c("TP53", "INS", "EGFR", "MET", NA)
#'   mget2(sym, org.Hs.egSYMBOL2EG)
#'   mget2(sym, org.Hs.egSYMBOL2EG, na.rm=FALSE)
#' }
#' 
mget2 <- function(x, y, sort=TRUE, first=TRUE, na.rm=TRUE) {
	idx <- !is.na(x)
	length(idx)>0 || return(NA)
	
	res <- rep(NA, length(x))
	res[idx] <- mget(x[idx], y, ifnotfound=NA)
	if( sort )  res <- lapply(res, sort)
	if( first ) res <- lapply(res, "[", 1)
	if( na.rm ) res <- res[!is.na(res)]
	
	res <- unlist(res)
	res	
}
# 2011-12-30
# - allow NA's in x
# 2012-03-06
# - removed dependence on pwbc::na.rm, for improved portability
