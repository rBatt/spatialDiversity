#' ---
#' title: "Spatial clustering of richness change in coastal North America"
#' author: "Ryan Batt"
#' date: "2017-02-03"
#' abstract: |
#'        North American coastal regions show long-term changes in species abundance and distribution, and changes in community diversity. What are the spatial patterns underlying these dynamics? Specifically, do changes in biodiversity occur homogenously or randomly across space, or are these proceses spatially patterned?
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     fig_caption: true
#'     theme: "readable"
#'     template: default
#'   pdf_document:
#'     toc: true
#'     toc_depth: 4
#'     template: latex-ryan.template
#'     fig_caption: true
#' geometry: margin=1.0in
#' lineno: true
#' lineSpacing: false
#' titlesec: true
#' documentclass: article
#' placeins: true
#' ---

#+ deleted-pandoc-headers, include=FALSE, echo=FALSE
# #'      pandoc_args: [
# #'      "--chapters"
# #'      ]


#+ setup, include=FALSE, echo=FALSE
# =================
# = Load Packages =
# =================
# #' date: "2016-08-10"
# #' title: "manuscript_results.R"
library('trawlDiversity')
library('spatialDiversity')
library('rbLib')
library("data.table")
library('maps')
library('raster')
library('spatstat')
library("spdep")

# Report
library(knitr)
library(rmarkdown)
# library(xtable)
# library(kfigr) # devtools::install_github("github mkoohafkan/kfigr") # ?
# library(stargazer)
# library(texreg)


# ================
# = Report Setup =
# ================
doc_type <- c("html", "pdf")[1]
table_type <- c("html"="html", "pdf"="latex")[doc_type]
options("digits"=3) # rounding output to 4 in kable() (non-regression tables)
o_f <- paste(doc_type, "document", sep="_")

# problem with pdflatex in El Capitan? It might be directories. Check http://pages.uoregon.edu/koch/FixLink.pkg

# render!
# rmarkdown::render(
# 	"~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript/reportSD.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript',
# 	clean = TRUE
# )


Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))

opts_chunk$set(
	fig.path = 'manuscript_report/', 
	cache.path='manuscript_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=F,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


# setwd("~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript")
source("../manuscript/manuscript_figures_functions.R")
# source("../manuscript/fig_tbl_number.R")
# eval(fig_tbl_number())


# ============
# = Richness =
# ============
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Data Setup
#+ dataSetup, echo=TRUE
mapDat <- spatialDiversity::mapDat
ureg <- mapDat[,unique(reg)]
#' 
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Spatial Clustering of Colonization and Extinction
#' ###Heat Maps
#' ####Figure 1. Richness map
#+ rich-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 1.** Maps of long-term averages of richness at each site for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Richness values were smoothed using a Gaussian kernel smoother. The smoothed richness value is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="richness", main="Average Richness")
#' 
# #' ####Figure 2. Colonization map
# #+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2.** Maps of long-term averages of colonizations per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of colonization rate were smoothed using a Gaussian kernel smoother. The smoothed colonization rate is indicated by the color bars in each panel; colors are scaled independently for each region."
# ceRate_map(ce="colonization")
# #'
# #' ####Figure 2b. Unique Colonization map
# #+ ucol-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2b.** Maps of the number of unique species with regional colonizations involving each site."
# ceRate_map(ce="uCol")
#' 
#' ####Figure 2c. Total Colonization map
#+ totcol-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2c.** Maps of total colonizations per site."
ceRate_map(ce="totCol", main="Total Colonizations")
#'   
# #' ####Figure 3. Extinction map
# #+ ext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 3.** Maps of long-term averages of extinctions per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of extinction rate were smoothed using a Gaussian kernel smoother. The smoothed extinction rate is indicated by the color bars in each panel; colors are scaled independently for each region."
# ceRate_map(ce="extinction")
# #'
# #' ####Figure 3b. Unique Extinction map
# #+ uext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 3b.** Maps of the number of unique species with regional colonizations involving each site."
# ceRate_map(ce="uExt")
#'   
#' ####Figure 3c. Total Extinction map
#+ totext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 3c.** Maps of the total number of regional extinctions involving each site."
ceRate_map(ce="totExt", main="Total Extinctions")
#' 
#' Hotspots can be seen in most regions. Newfoundland also has high values around its edge (as opposed to interior), it seems. NEUS and Gmex show very strong hotspots, and other locations tend to be much much lower. Other regions show more of a continuum.  
#'     
#+ col-ext-intensities, echo=TRUE,  cache=FALSE
sppp <- function(...){spatstat::Smooth(spatstat::ppp(...), hmax=1)}
map_smooth <- function(X, val=c("n_spp_col_weighted","n_spp_ext_weighted","avgRich","uCol","uExt","totCol","totExt")){
	val <- match.arg(val)
	r <- X[,unique(reg)]
	sppp(x=X[,lon], y=X[,lat], marks=X[,get(val)], window=mapOwin[[r]])
}
rel_col_ext_rate <- mapDat[,j={
	map_smooth_rich <- map_smooth(.SD, "avgRich")
	mark_range_rich <- range(map_smooth_rich, na.rm=TRUE)*10
	
	map_smooth_col <- map_smooth(.SD, "n_spp_col_weighted")
	mark_range_col <- range(map_smooth_col, na.rm=TRUE)*10
	
	map_smooth_ext <- map_smooth(.SD, "n_spp_ext_weighted")
	mark_range_ext <- range(map_smooth_ext, na.rm=TRUE)*10
	
	map_smooth_uCol <- map_smooth(.SD, "uCol")
	mark_range_uCol <- range(map_smooth_uCol, na.rm=TRUE)*10
	
	map_smooth_uExt <- map_smooth(.SD, "uExt")
	mark_range_uExt <- range(map_smooth_uExt, na.rm=TRUE)*10
	
	map_smooth_totCol <- map_smooth(.SD, "totCol")
	mark_range_totCol <- range(map_smooth_totCol, na.rm=TRUE)*10
	
	map_smooth_totExt <- map_smooth(.SD, "totExt")
	mark_range_totExt <- range(map_smooth_totExt, na.rm=TRUE)*10
	
	ol <- list(
		minval_rich=mark_range_rich[1], maxval_rich=mark_range_rich[2], 
		max_o_min_rich=do.call("/",as.list(rev(mark_range_rich))),
		minval_col=mark_range_col[1], maxval_col=mark_range_col[2], 
		max_o_min_col=do.call("/",as.list(rev(mark_range_col))),
		minval_ext=mark_range_ext[1], maxval_ext=mark_range_ext[2], 
		max_o_min_ext=do.call("/",as.list(rev(mark_range_ext))),
		minval_uCol=mark_range_uCol[1], maxval_uCol=mark_range_uCol[2], 
		max_o_min_uCol=do.call("/",as.list(rev(mark_range_uCol))),
		minval_uExt=mark_range_uExt[1], maxval_uExt=mark_range_uExt[2], 
		max_o_min_uExt=do.call("/",as.list(rev(mark_range_uExt))),
		minval_totCol=mark_range_totCol[1], maxval_totCol=mark_range_totCol[2], 
		max_o_min_totCol=do.call("/",as.list(rev(mark_range_totCol))),
		minval_totExt=mark_range_totExt[1], maxval_totExt=mark_range_totExt[2], 
		max_o_min_totExt=do.call("/",as.list(rev(mark_range_totExt)))
	)
	lapply(ol, function(x)if(is.numeric(x)){signif(x,3)}else{x})
},by=c("reg"), .SDcols=names(mapDat)]
# #+ col-ext-intensities-table, echo=FALSE
# kable(
# 	rbind(rel_col_ext_rate, rel_col_ext_rate[,lapply(.SD, median)][,reg:="MEDIAN"]),
# 	caption="The colonization and extinction intensity range and max/min ratio, and median among regions. Useful for assessing how big of a difference there is between red and blue for each region."
# )
#'   
#' ###Neighborhoods and Local Moran's I
#' ####Figure 4. Richness neighborhood
#+ rich-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 4.** Connectivity and local spatial autocorrelation of richness in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
nb_moranI(ce="richness")
#' 
# #' ####Figure 5. Colonization neighborhood
# #+ col-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 5.** Connectivity and local spatial autocorrelation of colonization events in each region."
# nb_moranI(ce="colonization")
# #'
# #' ####Figure 5b. Unique Colonization neighborhood
# #+ ucol-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 5b.** Connectivity and local spatial autocorrelation of colonization events (each species only counted once per stratum) in each region."
# nb_moranI(ce="uCol")
#' 
#' ####Figure 5c. Total Colonization neighborhood
#+ totcol-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 5c.** Connectivity and local spatial autocorrelation of colonization events (each species possibly counted more than once per stratum) in each region."
nb_moranI(ce="totCol")
#'   
# #' ####Figure 6. Extinction neighborhood
# #+ ext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 6.** Connectivity and local spatial autocorrelation of extinction events in each region."
# nb_moranI(ce="extinction")
# #'
# #' ####Figure 6b. Unique Extinction neighborhood
# #+ uext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 6b.** Connectivity and local spatial autocorrelation of extinction events (each species only counted once per stratum) in each region."
# nb_moranI(ce="uExt")
#'   
#' ####Figure 6c. Total Extinction neighborhood
#+ totext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 6c.** Connectivity and local spatial autocorrelation of extinction events (each species possibly counted more than once per stratum) in each region."
nb_moranI(ce="totExt")
#'   
#' \FloatBarrier  
#'   
#' ***  
#' 
#' ###Scatterplots involving colonization, richness, extinction
#' ####Figure 7. Richness vs Depth
#+ richVdepth, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 7.** Long-term average of per-site species richness vs the depth (m) of the site. Fitted line is a regression of richness ~ depth + depth^2."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(depth, avgRich);
		mtext(ureg[r],side=3,line=0.5,font=2);
		lines(sort(depth),predict(lm(avgRich~depth+I(depth^2),data=.SD[order(depth)])))
	}]
}

# par(mfrow=c(3,3));mapDat[,j={plot(depth, uCol);mtext(reg,side=3,line=0.5,font=2)},by='reg'] # no relationship between depth and number of species that had a colonization event involving the stratum, nor, as shown by switch uCol to uExt, between depth and extinction
#' 
#' ####Figure 8. Total Colonization vs Richness
#+ colVrich, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 8.** The total number of regional colonizations that involved a site vs long-term average of the site's richness."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(avgRich, totCol)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}
#' 
#' ####Figure 9.  Total Extinction vs Richness
#+ extVrich, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 9.** The total number of regional extinctions that involved a site vs long-term average of the site's richness."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(avgRich, totExt)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}

#' 
# #' ####Figure 10. Unique Colonization vs Unique Extinction
# #+ u-colVext, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 10.** Numbers of colonizations versus extinctions at each site. Species were only counted once per site (even if that species colonized or went extinct from the region multiple times, each time involving the site). Colors and shapes indicate the metrics with significant local spatial autocorrelation at the site: Blue +’s are colonization only, red x’s are extinction only, purple diamonds are both colonization and extinction, and black circles are neither colonization nor extinction. Gray circles were overlaid on sites with significant clustering in local species richness."
# eval(figure_setup())
# par(mfrow=c(3,3), mar=c(2.15,2.15,1.15,0.5), cex=1, mgp=c(1,0.25,0), tcl=-0.15, ps=10)
# for(r in 1:length(ureg)){
# 	mapDat[reg==ureg[r],j={
# 		sigColInd <- lI_pvalue_uCol<0.05
# 		sigExtInd <- lI_pvalue_uExt<0.05
# 		muCol <- mean(uCol)
# 		muExt <- mean(uExt)
# 		hotspotIndCol <- sigColInd #& (totCol > muCol)
# 		hotspotIndExt <- sigExtInd #& (totExt > muExt)
# 		both <- hotspotIndExt&hotspotIndCol
# 		neither <- !hotspotIndExt&!hotspotIndCol
#
# 		cols <- vector("character", length(hotspotIndCol))
# 		cols[hotspotIndCol] <- "blue"
# 		cols[hotspotIndExt] <- "red"
# 		cols[both] <- "purple"
# 		cols[neither] <- "black"
#
# 		pchs <- vector("integer", length(hotspotIndCol))
# 		pchs[hotspotIndCol] <- 3
# 		pchs[hotspotIndExt] <- 4
# 		pchs[both] <- 5
# 		pchs[neither] <- 1
#
# 		# print(all(hotspotIndCol | hotspotIndExt | both | neither))
#
# 		# pchs <- rep(21, length(hotspotIndExt))
# 		# pchs[!neither] <- 19
#
# 		plot(uExt, uCol, col=cols, pch=pchs, cex=1.2)
# 		mtext(pretty_reg[ureg[r]],side=3,line=0.01,font=2)
#
# 		sigRichInd <- lI_pvalue_rich<0.05
# 		hotspotIndRich <- sigRichInd
# 		points(uExt[hotspotIndRich], uCol[hotspotIndRich], col='gray', pch=20, cex=0.7)
# 	}]
# }
#' ####Figure 10b. Total Colonization vs Total Extinction
#+ tot-colVext, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 10b.** Total numbers of colonizations versus extinctions at each site. Colors and shapes indicate the metrics with significant local spatial autocorrelation at the site: Blue +’s are colonization only, red x’s are extinction only, purple diamonds are both colonization and extinction, and black circles are neither colonization nor extinction. Gray circles were overlaid on sites with significant clustering in local species richness."
eval(figure_setup())
par(mfrow=c(3,3), mar=c(2.15,2.15,1.15,0.5), cex=1, mgp=c(1,0.25,0), tcl=-0.15, ps=10)
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		sigColInd <- lI_pvalue_totCol<0.05
		sigExtInd <- lI_pvalue_totExt<0.05
		
		muCol <- mean(totCol)
		muExt <- mean(totExt)
		hotspotIndCol <- sigColInd #& (totCol > muCol)
		hotspotIndExt <- sigExtInd #& (totExt > muExt)
		
		both <- hotspotIndExt&hotspotIndCol
		neither <- !hotspotIndExt&!hotspotIndCol
		
		# print(all(hotspotIndCol | hotspotIndExt | both | neither))
		
		cols <- vector("character", length(hotspotIndCol))
		cols[hotspotIndCol] <- "blue"
		cols[hotspotIndExt] <- "red"
		cols[both] <- "purple"
		cols[neither] <- "black"
		
		pchs <- vector("integer", length(hotspotIndCol))
		pchs[hotspotIndCol] <- 3
		pchs[hotspotIndExt] <- 4
		pchs[both] <- 5
		pchs[neither] <- 1
		
		# pchs <- rep(21, length(hotspotIndExt))
		# pchs[!neither] <- 19
		
		plot(totExt, totCol, col=cols, pch=pchs, cex=1.2)
		mtext(pretty_reg[ureg[r]],side=3,line=0.01,font=2)
		
		sigRichInd <- lI_pvalue_rich<0.05
		hotspotIndRich <- sigRichInd
		points(totExt[hotspotIndRich], totCol[hotspotIndRich], col='gray', pch=20, cex=0.7)
	}]
}
#' 
#' ####Figure 11. Total Colonization vs Unique Colonization
#+ totcolVucol, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 11.** The total number of colonization events at each site vs the number of species that ever had a colonization event involving the site."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(uCol, totCol)
		abline(a=0, b=1)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}
#' 
#' ####Figure 12. Total Extinction vs Unique Extinction
#+ totextVuext, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure 12.** The total number of extinction events at each site vs the number of species that ever had an extinction event involving the site."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(uExt, totExt)
		abline(a=0, b=1)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}
#'   
#' \FloatBarrier  
#'   
#' ***  
#' 
#+ tbl-fracColExtInHotspot
fracHotspot <- mapDat[,j={
	sigColInd <- lI_pvalue_totCol<0.05
	sigExtInd <- lI_pvalue_totExt<0.05
	muCol <- mean(totCol)
	muExt <- mean(totExt)
	hotspotIndCol <- sigColInd & (totCol > muCol)
	hotspotIndExt <- sigExtInd & (totExt > muExt)
	
	list(
		hotspotFracCol=sum(totCol[hotspotIndCol])/sum(totCol), 
		hotspotFracColNSPots=sum(hotspotIndCol)/length(stratum),
		hotspotFracExt=sum(totExt[hotspotIndExt])/sum(totExt), 
		hotspotFracExtNSPots=sum(hotspotIndExt)/length(stratum)
	)
},by=c("reg")]
kable(fracHotspot)
#'   
#' \FloatBarrier  
#'   
#' ***  
#' 