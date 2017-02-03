#' Figure Setup for Trawl
#' 
#' Generates region names and colors for regional trawl figures
#' 
#' @param ... not used
#' 
#' @value Returns a quoted expression that can be evaluated with eval()
#' 
#' @examples
#' local({
#' 	print(ls())
#' 	eval(figure_setup())
#' 	sapply(mget(ls()), print)
#' 	par(mfrow=c(3,3))
#' 	for(i in 1:9){
#' 		plot(1, pch=19, col=pretty_col[i], cex=2)
#' 		mtext(pretty_reg[i], side=3, line=-2, font=2)
#' 	}
#' })
figure_setup <- function(...){
	bquote({
		regs <- trawlDiversity::comm_master[,una(reg)] #sapply(p, function(x)x$processed[,una(reg)])
		pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
		pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
		names(pretty_col) <- names(pretty_reg)
	})
}

#' Panel Label
#' 
#' Function that labels panels
#' 
#' @param ... other arguments, currently not used
#' 
#' @return Returns nothing, but adds letter labels to panels
panLab <- function(...){
	pm <- par("mfg")
	nmat <- matrix(1:(prod(pm[3:4])), nr=pm[3], nc=pm[4])
	panelLabel <- LETTERS[nmat[pm[1], pm[2]]]
	return(panelLabel)
}

