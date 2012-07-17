barplot.topTable.F <- function(tt, rma, coef.colour=NULL, rma.colour=NULL, number=10) {
	par(las=2, mar=c(12,4,2,2))
	layout(matrix(c(1,2), ncol=2), widths=c(2,1))
	
	# where are the coefs that are hidden in the tt...
	start.idx <- 2
	end.idx <- min(unlist(mgrep(c("F", "AveExpr"), colnames(tt)))) - 1
	coefs <- tt[,c(start.idx:end.idx)]
	
	ids <- tt$ID[1:number]
	for(id in ids) {
		barplot(as.numeric(rma[id,]), main=id, col=rma.colour, names.arg=colnames(rma))
		barplot(as.numeric(coefs[tt$ID == id,]), main="Coefficients", col=coef.colour, names.arg=colnames(coefs))
		abline(h=0)
	}
}

barplot.topTable.contrasts.F <- function(tt, preFit, rma, coef.colour=NULL, rma.colour=NULL, preFit.col=NULL, number=10) {
	opar <- par(no.readonly=TRUE)
	par(las=2, mar=c(12,4,2,2), oma=c(2,0,0,0))
	layout(matrix(c(1,2,3), nrow=1), widths=c(2,1,1))
	
	# where are the coefs that are hidden in the tt...
	start.idx <- 2
	end.idx <- min(unlist(mgrep(c("F", "AveExpr"), colnames(tt)))) - 1
	coefs <- tt[,c(start.idx:end.idx)]
	
	ids <- tt$ID[1:number]
	for(id in ids) {
		barplot(as.numeric(rma[id,]), main=id, col=rma.colour, names.arg=colnames(rma))
		barplot(as.numeric(preFit$coefficients[id, ]), 
				main="Coefficients prior to fitContrasts", col=preFit.col,
				names.arg=colnames(preFit$coefficients))
		abline(h=0)
		barplot(as.numeric(coefs[tt$ID == id,]), main="Coefficients", col=coef.colour, names.arg=colnames(coefs))
		abline(h=0)
		
		label <- paste("ProbeSetID:", id, "F test: P=",prettyNum(tt$P.Value[tt$ID==id]),
						" FDR=", prettyNum(tt$adj.P.Val[tt$ID==id]))
		cat(label, "\n")
		mtext(label, side=1, outer=TRUE, line=1, las=1)
	}
	par(opar)
}

