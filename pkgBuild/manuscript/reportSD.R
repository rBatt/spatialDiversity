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
#' #Spatial Clustering of Colonization and Extinction
#' ##Figure 1. Richness map
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
#' ##Figure 2. Total Colonization map
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
#' ##Figure 3. Total Extinction map
#+ totext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 3c.** Maps of the total number of regional extinctions involving each site."
ceRate_map(ce="totExt", main="Total Extinctions")
#' 
#' Hotspots can be seen in most regions. Newfoundland also has high values around its edge (as opposed to interior), it seems. NEUS and Gmex show very strong hotspots, and other locations tend to be much much lower. Other regions show more of a continuum.  
#'     
#' ##Table. Relative intensities of col/ ext in maps
#+ col-ext-intensities, echo=TRUE,  cache=FALSE
sppp <- function(...){spatstat::Smooth(spatstat::ppp(...), hmax=1)}
map_smooth <- function(X, val=c("n_spp_col_weighted","n_spp_ext_weighted","avgRich","uCol","uExt","totCol","totExt")){
	val <- match.arg(val)
	r <- X[,unique(reg)]
	sppp(x=X[,lon], y=X[,lat], marks=X[,get(val)], window=mapOwin[[r]])
}
rel_col_ext_rate <- mapDat[,j={
	map_smooth_rich <- map_smooth(.SD, "avgRich")
	mark_range_rich <- range(map_smooth_rich, na.rm=TRUE)#*10
	
	# map_smooth_col <- map_smooth(.SD, "n_spp_col_weighted")
	# mark_range_col <- range(map_smooth_col, na.rm=TRUE)*10
	#
	# map_smooth_ext <- map_smooth(.SD, "n_spp_ext_weighted")
	# mark_range_ext <- range(map_smooth_ext, na.rm=TRUE)*10
	#
	# map_smooth_uCol <- map_smooth(.SD, "uCol")
	# mark_range_uCol <- range(map_smooth_uCol, na.rm=TRUE)*10
	#
	# map_smooth_uExt <- map_smooth(.SD, "uExt")
	# mark_range_uExt <- range(map_smooth_uExt, na.rm=TRUE)*10
	
	map_smooth_totCol <- map_smooth(.SD, "totCol")
	mark_range_totCol <- range(map_smooth_totCol, na.rm=TRUE)*10
	
	map_smooth_totExt <- map_smooth(.SD, "totExt")
	mark_range_totExt <- range(map_smooth_totExt, na.rm=TRUE)*10
	
	ol <- list(
		minval_rich=mark_range_rich[1], maxval_rich=mark_range_rich[2], 
		max_o_min_rich=do.call("/",as.list(rev(mark_range_rich))),
		# minval_col=mark_range_col[1], maxval_col=mark_range_col[2],
		# max_o_min_col=do.call("/",as.list(rev(mark_range_col))),
		# minval_ext=mark_range_ext[1], maxval_ext=mark_range_ext[2],
		# max_o_min_ext=do.call("/",as.list(rev(mark_range_ext))),
		# minval_uCol=mark_range_uCol[1], maxval_uCol=mark_range_uCol[2],
		# max_o_min_uCol=do.call("/",as.list(rev(mark_range_uCol))),
		# minval_uExt=mark_range_uExt[1], maxval_uExt=mark_range_uExt[2],
		# max_o_min_uExt=do.call("/",as.list(rev(mark_range_uExt))),
		minval_totCol=mark_range_totCol[1], maxval_totCol=mark_range_totCol[2], 
		max_o_min_totCol=do.call("/",as.list(rev(mark_range_totCol))),
		minval_totExt=mark_range_totExt[1], maxval_totExt=mark_range_totExt[2], 
		max_o_min_totExt=do.call("/",as.list(rev(mark_range_totExt)))
	)
	lapply(ol, function(x)if(is.numeric(x)){signif(x,3)}else{x})
},by=c("reg"), .SDcols=names(mapDat)]
#+ col-ext-intensities-table, echo=FALSE
kable(
	rbind(rel_col_ext_rate, rel_col_ext_rate[,lapply(.SD, median)][,reg:="MEDIAN"]),
	caption="The colonization and extinction intensity range and max/min ratio, and median among regions. Useful for assessing how big of a difference there is between red and blue for each region."
)
#'   
#' #Neighborhoods and Local Moran's I
#' ##Figure S1. Richness neighborhood
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
#' ##Figure S2. Total Colonization neighborhood
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
#' ##Figure S3. Total Extinction neighborhood
#+ totext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S3.** Connectivity and local spatial autocorrelation of extinction events (each species possibly counted more than once per stratum) in each region."
nb_moranI(ce="totExt")
#'   
#' \FloatBarrier  
#'   
#' ***  
#' 
#' #Scatterplots involving colonization, richness, extinction
#' ##Figure 4 & Table: Total Colonization vs Total Extinction
#+ tot-colVext, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure Col v Ext** Total numbers of colonizations versus extinctions at each site. Colors and shapes indicate the metrics with significant local spatial autocorrelation at the site: Blue +’s are colonization only, red x’s are extinction only, purple diamonds are both colonization and extinction, and black circles are neither colonization nor extinction. Gray circles were overlaid on sites with significant clustering in local species richness."
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
		
		plot(totExt[order(pchs)], totCol[order(pchs)], 
			pch=pchs[order(pchs)], col=cols[order(pchs)], cex=1.2,
			xlab="Extinctions", ylab="Colonizations"
		)
		abline(v=muExt, lty='dashed')
		abline(h=muCol, lty='dashed')
		mtext(pretty_reg[ureg[r]],side=3,line=0.01,font=2)
		
		sigRichInd <- lI_pvalue_rich<0.05
		hotspotIndRich <- sigRichInd
		points(totExt[hotspotIndRich], totCol[hotspotIndRich], col='gray', pch=20, cex=0.7)
	}]
}

#+ col-ext-regression-table
colExt_regs <- list()
ureg <- mapDat[,unique(reg)]
for(r in 1:length(ureg)){
	colExt_regs[[r]] <- summary(lm(totCol~totExt, data=mapDat[reg==ureg[r]]))
}
names(colExt_regs) <- ureg
mod_sum_funct <- function(x){
	as.list(c(x$coefficients[2,c("Estimate","Pr(>|t|)")],"R2"=x$r.squared))
}
colExt_mod_sum <- lapply(colExt_regs, mod_sum_funct)
colExt_mod_sum <- cbind(reg=names(colExt_regs), rbindlist(colExt_mod_sum))
colExt_mod_sum[,c("pAdjusted"):=p.adjust(get("Pr(>|t|)"), method="BH")]
kable(colExt_mod_sum, caption="Table. Statistics for totCol~totExt linear regression.")

#' ##Figure: Richness vs Depth
#+ richVdepth, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure Rich v Depth** Long-term average of per-site species richness vs the depth (m) of the site. Fitted line is a regression of richness ~ depth + depth^2."
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
#' ##Figure & Table: Total Colonization vs Richness
#+ colVrich, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure Col v Rich** The total number of regional colonizations that involved a site vs long-term average of the site's richness."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(avgRich, totCol)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}

#+ colVrich-regression, results="asis"
colRich_regs <- list()
ureg <- mapDat[,unique(reg)]
for(r in 1:length(ureg)){
	colRich_regs[[r]] <- summary(lm(totCol~avgRich, data=mapDat[reg==ureg[r]]))
}
names(colRich_regs) <- ureg
mod_sum_funct <- function(x){
	as.list(c(x$coefficients[2,c("Estimate","Pr(>|t|)")],"R2"=x$r.squared))
}
colRich_mod_sum <- lapply(colRich_regs, mod_sum_funct)
colRich_mod_sum <- cbind(reg=names(colRich_regs), rbindlist(colRich_mod_sum))
colRich_mod_sum[,c("pAdjusted"):=p.adjust(get("Pr(>|t|)"), method="BH")]
kable(colRich_mod_sum, caption="Table. Statistics for totCol~avgRich linear regression.")

#' 
#' ##Figure & Table: Total Extinction vs Richness
#+ extVrich, echo=TRUE, fig.width=7, fig.height=7, fig.cap="**Figure Ext v Rich** The total number of regional extinctions that involved a site vs long-term average of the site's richness."
par(mfrow=c(3,3))
for(r in 1:length(ureg)){
	mapDat[reg==ureg[r],j={
		plot(avgRich, totExt)
		mtext(ureg[r],side=3,line=0.5,font=2)
	}]
}

#+ extVrich-regression
extRich_regs <- list()
ureg <- mapDat[,unique(reg)]
for(r in 1:length(ureg)){
	extRich_regs[[r]] <- summary(lm(totExt~avgRich, data=mapDat[reg==ureg[r]]))
}
names(extRich_regs) <- ureg
mod_sum_funct <- function(x){
	as.list(c(x$coefficients[2,c("Estimate","Pr(>|t|)")],"R2"=x$r.squared))
}
extRich_mod_sum <- lapply(extRich_regs, mod_sum_funct)
extRich_mod_sum <- cbind(reg=names(extRich_regs), rbindlist(extRich_mod_sum))
extRich_mod_sum[,c("pAdjusted"):=p.adjust(get("Pr(>|t|)"), method="BH")]
kable(extRich_mod_sum, caption="Table. Statistics for totExt~avgRich linear regression.")

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


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Compare Total and Unique Metrics
#' ##Figure: Total Colonization vs Unique Colonization
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
#' ##Figure: Total Extinction vs Unique Extinction
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
nLSA <- mapDat[,j={
	sigRichInd <- lI_pvalue_rich<0.05
	sigColInd <- lI_pvalue_totCol<0.05
	sigExtInd <- lI_pvalue_totExt<0.05
	muCol <- mean(totCol)
	muExt <- mean(totExt)
	hotspotIndCol <- sigColInd & (totCol > muCol)
	hotspotIndExt <- sigExtInd & (totExt > muExt)
	
	list(
		nSites = length(unique(stratum)),
		
		nRichLSA = sum(sigRichInd),
		nRichHotspot = sum(sigRichInd & (avgRich > mean(avgRich))),
		nRichColdspot = sum(sigRichInd & (avgRich < mean(avgRich))),
		
		nColLSA = sum(sigColInd),
		nColHotspot = sum(hotspotIndCol),
		nColColdspot = sum(sigColInd & (totCol < muCol)),
		
		nExtLSA = sum(sigExtInd),
		nExtHotspot = sum(hotspotIndExt),
		nExtExtdspot = sum(sigExtInd & (totExt < muExt))
	)
},by=c("reg")]
kable(nLSA)
#'   
#' \FloatBarrier  
#'   
#' ***  
#' 
#' ##Exploring 'Endemism' Patterns
#' One of the hypotheses was that sites with higher richness might have higher local colonizations and extinctions. If the richness at the site is also endemic -- meaning the species occurs at no other sites in the region -- then any local colonization or extinction of that endemic species is also a regional colonization or extinction. Therefore, in addition to knowing whether or not richness is correlated with colonization or extinction, we should also determine whether there is a degree of endemism at any of these sites. Otherwise, the simple proposed mechanism explaining relationships between richness and col/ext does not apply.
#+ endemism-calc, echo=TRUE
#' For each site in a region, how many of that site's occupants tend to occupy other sites at the same time?
endo <- trawlDiversity::data_all[reg!='wcann' & K==1,j={
	
	totStrat <- length(unique(stratum))
	list(stratum=stratum, totStrat=totStrat)
	
},by=c('reg','year','spp')]
print(endo[,sum(totStrat==1),by=c('reg','spp')])


#' Number of species that are only present in 1 site per year per region: `r endo[,max(totStrat),by=c('reg','spp')][,sum(V1==1)]`  
#' A list of the species that appear in only 1 site at a time:  
kable(endo[,max(totStrat),by=c('reg','spp')][V1==1])
#' Most of these examples are from the NEUS and GMEX. The one from SA is the Gaftopsail Sea Catfish. It's not super rare at all. In the full data set it was found in multiple strata in the same year, just not necessarily for K==1.  
#'   
#' Here's some examples of what I looked at from the NEUS  
sppImg("Syacium papillosum") 
data_all[spp=="Syacium papillosum", plot(lon, lat, col=as.factor(reg))]; map(add=TRUE)
#' This species is only ever found in 1 site at a time in NEUS. According to fish base (http://www.aquamaps.org/receive.php?type_of_map=regular), it should be fairly rare north of Capte Hatteras
#'   
sppImg("Decapterus punctatus")
data_all[spp=="Decapterus punctatus", plot(lon, lat, col=as.factor(reg))]; map(add=TRUE) 
#' Example of a species that only ever appeared in one site at a time (in a region). However, this strikes me as odd because the species is apparently more widely distributed throughout the NEUS than this map shows, according to fishbase (http://www.aquamaps.org/receive.php?type_of_map=regular). It's also a pelagic, so initially I just thought this was a sampling issue.
#'   
#' Now to examine the question from a site perspective.
e1 <- endo[,list(avg_totStrat=mean(totStrat)),by=c("reg","year","stratum")]
e2 <- e1[,list(avg_totStrat=mean(avg_totStrat)),by=c("reg","stratum")]
par(mfrow=c(3,3))
merge(mapDat, e2)[,j={hist(avg_totStrat, main=reg[1]);NULL}, by='reg']


#' #Beta Diversity Clustering of Sites w/ Local AC
#' ##Geographic Distance and Beta Diversity Functions
#+ betaD-AC-Cluster-functions
betaDist <- function(Y){
	ade4::dist.binary(Y, method=1)^2
}

geoDist <- function(x, y){
	if(!is.null(nrow(x))){
		x <- as.matrix(x)
	}
	if(!is.null(nrow(y))){
		y <- as.matrix(y)
	}
	
	x180 <- x < -180
	x[x180] <- x[x180] + 360
	y180 <- y < -180
	y[y180] <- y[y180] + 360
	
	geosphere::distVincentyEllipsoid(x, y)/1E3
}

dendroMap <- function(Dat){
	ur <- Dat[,unique(reg)]
	mat <- structure(vector('list', length(ur)), .Names=ur)
	beta <- structure(vector('list', length(ur)), .Names=ur)
	geo <- structure(vector('list', length(ur)), .Names=ur)
	for(r in 1:length(ur)){
		tdt <- Dat[reg==ur[r],list(spp=unique(spp), pres=1),by=c('reg','stratum')]
		tmat <- reshape2::acast(tdt, stratum~spp, value.var='pres')
		tmat[is.na(tmat)] <- 0
		mat[[ur[r]]] <- tmat
	
		if(nrow(tmat) > 1){
			bD <- betaDist(tmat)
			hcbd <- hclust(bD); 
			clusts <- cutree(hcbd, h=0.7)
		
			col_opts <- viridis::viridis(length(unique(clusts)))
			nameCol_key <- data.table(stratum=hcbd$labels, ord=hcbd$order, short=LETTERS[1:nrow(tmat)], clusts=clusts, col=col_opts[clusts])
		
			strat_names <- nameCol_key[,stratum] #hcbd$labels[hcbd$order] #rownames(tmat)
			strat_names_short <- structure(nameCol_key[,short], .Names=nameCol_key[,stratum])
		
			clustCol <- nameCol_key[,unique(col)[unique(clusts[ord])]] #nameCol_key[,unique(col)] #viridis::viridis(length(unique(clusts)))
			stratCol <- nameCol_key[,col] #clustCol[clusts]
			names(stratCol) <- nameCol_key[,stratum] #strat_names
		
			hcbd$labels <- nameCol_key[,short]#[order(hcbd$order)] #strat_names_short[order(hcbd$order)]
		}else{
			if(nrow(tmat)>0){
				strat_names <- rownames(tmat)
				strat_names_short <- structure(LETTERS[1:nrow(tmat)], .Names=strat_names)
				clustCol <- viridis::viridis(1)
				stratCol <- clustCol
				names(stratCol) <- strat_names
			}
		}

		plot(mapOwin[[ur[r]]], main='')
		if(nrow(tmat)>0){
			ll <- strat2ll(strat_names)
			text(x=(ll$lon), y=(ll$lat), strat_names_short, col=stratCol[strat_names])
			# geo <- geoDist(ll$lon, ll$lat) # geoDist(c(ll$lon[1], ll$lat[1]), c(ll$lon[2], ll$lat[2]))
		}
	
		if(nrow(tmat) > 1){
			plot(as.dendrogram(hcbd), leaflab='perpendicular')
			if(sum(hcbd$height < 0.7)>1 & sum(hcbd$height > 0.7)>1){
				rect.hclust(hcbd, h=0.7, border=clustCol)
			}
		}else{
			plot(1,1, type='n', xlab='', ylab='', xaxt='n', yaxt='n', bty='l')
		}
	}
	invisible(NULL)
}

strat2ll <- function(stratum){
	ll_split <- strsplit(stratum, " ")
	Lon <- as.numeric(sapply(ll_split, function(x)x[1]))
	Lat <- as.numeric(sapply(ll_split, function(x)x[2]))
	list(lon=Lon, lat=Lat)
}



# ---- cluster sites ----
#' ##Data for Col/ Ext Clustering & Local Community Composition
#+ betaD-AC-Cluster-data
#' First get the sites that have significant local spatial AC  
colClustSite <- spatialDiversity::mapDat[lI_pvalue_totCol < 0.05 & totCol>mean(totCol), stratum, by='reg']
extClustSite <- spatialDiversity::mapDat[lI_pvalue_totExt < 0.05 & totExt>mean(totExt), stratum, by='reg']

#' Pull out the years and species for colonizations, extinctions  
colSppYear <- col_ext_dt[col==1, list(spp, year), by='reg']
extSppYear <- col_ext_dt[ext==1, list(spp, year), by='reg']

#' Cross reference the AC sites and the c/e years & spp to subset full data set  
#' Basically, had to use mapDat to get there "where", and had to use col_ext_dt to get the "who" and "when",  then I had to use data_all2 (same as data_all, except subset to first haul within stratum-year) to get wtcpue, btemp, depth, etc.
subCols <- expression(list(reg, year, spp, stratum, K, Kmax, lon, lat, wtcpue, btemp, depth))
colDat <- data_all2[colSppYear, on=c('reg','spp','year')][colClustSite, on=c("reg", "stratum"), eval(subCols)] #[!is.na(lon)]
extDat <- data_all2[extSppYear, on=c('reg','spp','year')][extClustSite, on=c("reg", "stratum"), eval(subCols)]#[!is.na(lon)]


# lay_mat0 <- as.matrix(raster::disaggregate(raster::raster(matrix(1:9, nrow=3)), fact=4))
# lay_mat <- lay_mat0
# for(r in 1:length(ur)){
# 	val <- 1 + 2*(r-1)
# 	ind <- lay_mat0==r
# 	lay_mat[ind] <- val
# 	nr <- sqrt(sum(ind))
# 	slm <- (lay_mat[ind])
# 	tm <- matrix(lay_mat[ind], nrow=nr)
# 	diag(tm) <- c(rep(val+1, nr/2), rep(val, nr/2))
# 	lay_mat[ind] <- tm
# }

#' ##Figure: Colonization AC & Beta Diversity Clustering
#+ betaD-AC-Cluster-colFig, fig.cap="**Exploratory Figure.** Sites are labeled with a letter if their rates of EXTINCTION have significant local spatial autocorrelation (AC). In the left-hand panels, their geographic locations are shown. In the right-hand panels are dendrograms that are drawn via hierarchical clustering using beta diversity as the distance metric. Thus, left-hand panels indicate clustering based on spatial AC of EXTINCTION rates, and the right-hand panels indicate clustering based on beta diversity. Color boxes are drawn around dendrograms clusters (though some do not show up, and I'm not entirely sure why)."
par(mfrow=c(9,2), mar=c(1,1,0.25,0.25), cex=1, ps=8)
dendroMap(colDat)

#' ##Figure: Extinction AC & Beta Diversity Clustering
#+ betaD-AC-Cluster-extFig, fig.cap="**Exploratory Figure.** Sites are labeled with a letter if their rates of COLONIZATION have significant local spatial autocorrelation (AC). In the left-hand panels, their geographic locations are shown. In the right-hand panels are dendrograms that are drawn via hierarchical clustering using beta diversity as the distance metric. Thus, left-hand panels indicate clustering based on spatial AC of COLONIZATION rates, and the right-hand panels indicate clustering based on beta diversity. Color boxes are drawn around dendrograms clusters (though some do not show up, and I'm not entirely sure why)."
par(mfrow=c(9,2), mar=c(1,1,0.25,0.25), cex=1, ps=8)
dendroMap(extDat)

# ur <- colDat[,unique(reg)]
# for(r in 1:length(ur)){
# 	tdt <- colDat[reg==ur[r],list(spp=unique(spp), pres=1),by=c('reg','stratum')]
# 	tmat <- reshape2::acast(tdt, stratum~spp, value.var='pres')
# 	tmat[is.na(tmat)] <- 0
#
# 	if(nrow(tmat) > 1){
# 		bD <- betaDist(tmat)
# 		hcbd <- hclust(bD);
# 		clusts <- cutree(hcbd, h=0.7)
#
# 		col_opts <- viridis::viridis(length(unique(clusts)))
# 		nameCol_key <- data.table(stratum=hcbd$labels, ord=hcbd$order, short=LETTERS[1:nrow(tmat)], clusts=clusts, col=col_opts[clusts])
#
# 		strat_names <- nameCol_key[,stratum] #hcbd$labels[hcbd$order] #rownames(tmat)
# 		strat_names_short <- structure(nameCol_key[,short], .Names=nameCol_key[,stratum])
#
# 		clustCol <- nameCol_key[,unique(col)[unique(clusts[ord])]] #nameCol_key[,unique(col)] #viridis::viridis(length(unique(clusts)))
# 		stratCol <- nameCol_key[,col] #clustCol[clusts]
# 		names(stratCol) <- nameCol_key[,stratum] #strat_names
#
# 		hcbd$labels <- nameCol_key[,short]#[order(hcbd$order)] #strat_names_short[order(hcbd$order)]
# 	}else{
# 		if(nrow(tmat)>0){
# 			strat_names <- rownames(tmat)
# 			strat_names_short <- structure(LETTERS[1:nrow(tmat)], .Names=strat_names)
# 			clustCol <- viridis::viridis(1)
# 			stratCol <- clustCol
# 			names(stratCol) <- strat_names
# 		}
# 	}
#
# 	plot(mapOwin[[ur[r]]], main='')
# 	if(nrow(tmat)>0){
# 		ll <- strat2ll(strat_names)
# 		text(x=(ll$lon), y=(ll$lat), strat_names_short, col=stratCol[strat_names])
# 	}
#
# 	if(nrow(tmat) > 1){
# 		plot(as.dendrogram(hcbd), leaflab='perpendicular')
# 		if(sum(hcbd$height < 0.7)>1 & sum(hcbd$height > 0.7)>1){
# 			rect.hclust(hcbd, h=0.7, border=clustCol)
# 		}
# 	}else{
# 		plot(1,1, type='n', xlab='', ylab='', xaxt='n', yaxt='n', bty='o')
# 	}
# }




