#' mget with multiple search keys
#' 
#' lookup keys in a heirarchy of mapping environments, returning a value
#' for each key in the earliest mapping environment possible.
#' 
#' For a set of microarray probe ID's, sometimes you want to determine the
#' best ID for that probe, starting with Gene Symbol, if known, falling back
#' to RefSeq ID, UniGene ID, Genbank ID, and finally, just the probe ID.
#' There is thus a heirarchy of annotation tables.
#' This function looks up the search keys in the first mapping environment. If no value
#' is found for some keys, then those keys are searched within the next environment
#' and so it continues until you've run out of mapping environments.
#'
#' @note This function uses the \code{Lkeys} of the mapping tables to chain them together.
#' What about an approach that uses the \code{Rkeys}, so you could get probe ID's from
#' a vector of heterogeneous ID types.
#' 
#' @param x a vector of search keys
#' @param \dots named search environments, passed to \code{\.merge.values}
#' @return a named vector of values, corresponding to each of the search keys
#' @author Mark Cowley, 2011-09-27
#' @export
#' @examples
#' # if( require(org.Hs.eg.db) ) {
#' #   ids <- c("TP53", "INS", "SST", "12CC4")
#' #   eg <- mget.multi(ids, revmap(org.Hs.egSYMBOL), org.Hs.egALIAS2EG)
#' #   eg
#' #   mget2(eg, org.Hs.egSYMBOL)
#' # }
#' if( require(org.Hs.eg.db) && require(illuminaHumanv4.db) ) {
#'   ids <- c("ILMN_1779356", "ILMN_1666966", "ILMN_1812824", "ILMN_3302350", "ILMN_1343059")
#'   mget.multi(ids, illuminaHumanv4SYMBOL, illuminaHumanv4REFSEQ, illuminaHumanv4ENSEMBL, illuminaHumanv4ACCNUM)
#' }
#'
mget.multi <- function(x, ...) {
	bestId <- .merge.values(...)
	res <- bestId[match(x, names(bestId))]
	
	res
}


#' @noRd
.merge.values <- function(...) {
	# TODO: use match.calls to get the map names.
	
	maps <- list(...)
	alleq(sapply(maps, Llength)) || stop("The Lkeys must match in all the provided mapping objects")
	
	keys <- Lkeys(maps[[1]])
	res <- rep(NA, length(keys))
	for(map in maps) {
		all(keys %in% Lkeys(map)) || stop(sprintf("Lkeys in map %s don't match the keys in the first map", map@objName))
		idx <- is.na(res)
		vals <- mget2(keys[idx], map, sort=FALSE, na.rm=FALSE)
		message(sprintf("Found %05d values using %s", sum(!is.na(vals)), map@objName))
		res[idx] <- vals
	}
	# if there's still NA's afterwards, use the Lkeys as the best ID
	idx <- is.na(res)
	res[idx] <- keys[idx]
	message(sprintf("Found %05d values using Lkeys of first map", sum(idx)))
	names(res) <- keys
	
	res
}