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
source("../../../trawlDiversity/pkgBuild/manuscript/manuscript_figures_functions.R")
# source("../manuscript/fig_tbl_number.R")
# eval(fig_tbl_number())


# ============
# = Richness =
# ============
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Spatial Clustering of Colonization and Extinction
#' ###Heat Maps
#' ####Figure 1. Richness map
#+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2.** Maps of long-term averages of richness at each site for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Richness values were smoothed using a Gaussian kernel smoother. The smoothed richness value is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="richness")
#' 
#' ####Figure 1. Colonization map
#+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2.** Maps of long-term averages of colonizations per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of colonization rate were smoothed using a Gaussian kernel smoother. The smoothed colonization rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="colonization")
#'   
#' ####Figure 2. Extinction map
#+ ext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S3.** Maps of long-term averages of extinctions per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of extinction rate were smoothed using a Gaussian kernel smoother. The smoothed extinction rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="extinction")
#' 
#' Hotspots can be seen in most regions. Newfoundland also has high values around its edge (as opposed to interior), it seems. NEUS and Gmex show very strong hotspots, and other locations tend to be much much lower. Other regions show more of a continuum.  
#'     
#+ col-ext-intensities, echo=TRUE,  cache=FALSE
sppp <- function(...){spatstat::Smooth(spatstat::ppp(...), hmax=1)}
map_smooth <- function(X, val=c("n_spp_col_weighted","n_spp_ext_weighted","avgRich")){
	val <- match.arg(val)
	r <- X[,unique(reg)]
	sppp(x=X[,lon], y=X[,lat], marks=X[,get(val)], window=mapOwin[[r]])
}
rel_col_ext_rate <- mapDat[,j={
	map_smooth_col <- map_smooth(.SD, "n_spp_col_weighted")
	mark_range_col <- range(map_smooth_col, na.rm=TRUE)*10
	
	map_smooth_ext <- map_smooth(.SD, "n_spp_ext_weighted")
	mark_range_ext <- range(map_smooth_ext, na.rm=TRUE)*10
	
	ol <- list(
		minval_col=mark_range_col[1], maxval_col=mark_range_col[2], 
		max_o_min_col=do.call("/",as.list(rev(mark_range_col))),
		minval_ext=mark_range_ext[1], maxval_ext=mark_range_ext[2], 
		max_o_min_ext=do.call("/",as.list(rev(mark_range_ext)))
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
#' ####Figure S4. Colonization neighborhood
#+ col-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S4.** Connectivity and local spatial autocorrelation of colonization events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
nb_moranI(ce="colonization")
#'   
#' ####Figure S5. Extinction neighborhood
#+ ext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S5.** Connectivity and local spatial autocorrelation of extinction events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
nb_moranI(ce="extinction")
#'   
#' \FloatBarrier  
#'   
#' ***  

