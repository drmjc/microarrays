# Compare two topTables via correlation of the t-stats
#
# Parameters:
#	tt1, tt2: two instances of topTables. They can have different numbers of probes, but some overlap is required.
# xlab, ylab, main: obvious
# column: which column to plot. default is "t". "logFC" would be another good choice.
#	...: additional arguments passed to plot.cor. see plot.cor
#
# Value:
#	creates a correlation plot. see plot.cor
#	invisibly returns the 2 toptables merged on a common ID. The unique suffixes on the colnames comes from xlab & ylab.
#
# See Also:
#	plot.cor
#
# Mark Cowley, 2010-01-18
# 2010-12-09: update: allow any column, not just t; change main; invisibly return the data
#


#' Compare two topTables via correlation of the t-stats
#' 
#' @param tt1 two instances of topTables. They can have different numbers of
#'   probes, but some overlap is required.
#' @param tt2 two instances of topTables. They can have different numbers of
#'   probes, but some overlap is required.
#' @param xlab obvious
#' @param ylab obvious
#' @param main obvious
#' @param column which column to plot. default is "t". "logFC" would be another
#'   good choice.
#' @param \dots additional arguments passed to plot.cor. see plot.cor
#' @return creates a correlation plot. see plot.cor invisibly returns the 2
#'   toptables merged on a common ID. The unique suffixes on the colnames comes
#'   from xlab & ylab.  See Also: plot.cor
#' @author Mark Cowley, 2010-01-18
#' @export
plot_cor_topTable <- function (tt1, tt2, xlab = "topTable 1", ylab = "topTable 2", main="topTable correlation", column="t", ...)  {
   ids <- intersect(tt1$ID, tt2$ID)
   tt1 <- tt1[match(ids, tt1$ID), ]
   tt2 <- tt2[match(ids, tt2$ID), ]
   plot.cor(tt1[,column], tt2[,column], main = main, xlab = xlab,
       ylab = ylab, ...)
   res <- merge(tt1, tt2, by="ID", suffixes=c(paste(".", xlab, sep=""), paste(".", ylab, sep="")))
   invisible(res)
}
