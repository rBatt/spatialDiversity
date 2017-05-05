
library(qgraph)
library(igraph)
library(trawlDiversity)
library(spatialDiversity)

setwd("~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript")
source("../manuscript/manuscript_figures_functions.R")

# =======================================
# = Core-Transient Cross-Region Network =
# =======================================
sp2 <- spp_master[,.SD[length(unique(reg))>1 & any(ce_categ!='neither') & any(ce_categ=='neither')],by=c('spp')]
relations <- sp2[,j={
	from_regs <- .SD[ce_categ=='neither', unique(reg)]
	to_regs <- .SD[ce_categ!='neither', unique(reg)]
	CJ(from_reg=from_regs, to_reg=to_regs)
	
},by=c("spp")]
# g <- graph_from_data_frame(relations[,list(from=from_reg, to=to_reg, spp=spp)], directed=TRUE)
g <- igraph::graph_from_edgelist(as.matrix(relations[,list(from=from_reg, to=to_reg)]), directed=TRUE)
adj_graph <- igraph::graph_from_adjacency_matrix(igraph::as_adj(g), mode='directed', weighted=TRUE)

# plot(g)
# plot(adj_graph, edge.width=(igraph::E(adj_graph)$weight)/3, edge.arrow.width=1.5, edge.arrow.size=1.25, edge.curved=TRUE)

# png("~/Desktop/trawlNetwork.png", res=200, units='in', width=4, height=4)
qgraph::qgraph(as_adj(g), colFactor=1, edge.width=1.5, labels=names(V(adj_graph)), threshold=1)
# dev.off()


# ==============================================
# = Colonization-Extinction Cross-Site Network =
# ==============================================
# Part 1 (spp_master)
# identify year of each colonization (per species-region)
# identify year of each extinction
# for each colonization year, if there is not a later year for extinction (or equal, if only present for 1 year in a row), remove it
# for each remaining colonization year, assign it the 'from' attribute
# for each colonization year that is 'from', identify an extinction year that is later, and assign it the 'to' attribute
# summarize result as a reg, spp, from_year, to_year
# 
# Part 2 (data_all)
# for each reg-spp-from_year combination, lookup the colonization sites, assign from_site
# for each reg-spp-to_year combination, lookup the extinction sites, assign to_site
# use CJ() to associate each from_site with each to_site for a given reg-spp-from_year-to_year combination
# summarize result as reg, spp, from_year, from_site, to_year, to_site

findSoonestAfter <- function(x, y){
	# @param x length 1, a year for which we want to find a subsequent or concurrent event 
	# @param y a vector of arbitrary length, containing years that may or may not occur after x
	#
	# @details
	# Finds the smallest value in y that is greater than or equal to x Finds indices of y whose values are greater than or equal to the value of x. Then
	stopifnot(length(x)==1)
	later <- y>=x
	if(any(later, na.rm=TRUE)){
		min(y[later], na.rm=TRUE)
	}else{
		NA
	}
	# would just do match(TRUE, y>=x), but this only returns earliest year (years are in y) if y is sorted
	# SLOWER:
	# fsa2 <- function(x, y){
# 		stopifnot(length(x)==1)
# 		y <- sort(y)
# 		later <- y>=x
# 		y[match(TRUE,later)]
# 	}
}


from_to_years <- trawlDiversity::spp_master[ce_categ=="both",j={
	col_years_all <- year[col==1]
	ext_years_all <- year[ext==1]
	
	to_year <- sapply(col_years_all, findSoonestAfter, y=ext_years_all)
	noNA <- !is.na(to_year)
	if(any(noNA)){
		from_year <- col_years_all[noNA]
		to_year <- to_year[noNA]
		list(from_year=from_year, to_year=to_year)
	}else{
		NULL # easier to do null than to have to complete.cases() later
		# list(from_year=NA_real_, to_year=NA_real_)
	}
},by=c("reg","spp")]


# This is slow ... see faster version below doign pt1, pt2 etc
# from_to_years[,j={
# 	from_sites <- data_all[.SD, on=c("reg","spp",year="from_year")][,stratum]
# 	to_sites <- data_all[.SD, on=c("reg","spp",year="to_year")][,stratum]
# 	combos <- CJ(from_site=from_sites, to_site=to_sites)
# 	data.table(from_year, to_year, combos)
# },by=c("reg","spp","from_year","to_year"),.SDcols=names(from_to_years)]

# faster than above
pt1 <- unique(trawlDiversity::data_all[,list(reg,spp,year,stratum)])[from_to_years,on=c("reg","spp",year="from_year")]
pt1 <- pt1[,list(reg,spp,from_year=year,to_year,from_site=stratum)]
pt2 <- unique(trawlDiversity::data_all[,list(reg,spp,year,stratum)])[from_to_years,on=c("reg","spp",year="to_year")]
pt2 <- pt2[,list(reg,spp,from_year,to_year=year,to_site=stratum)]
from_to <- merge(pt1, pt2, all=TRUE, allow.cartesian=TRUE)

stratum2ll <- function(x){
	lonlat_pat <- "^(-[0-9]{2,3}\\.[0-9]{2}) ([0-9]{2,3}\\.[0-9]{2}) [0-9]{1,4}$"
	data.table(site=x, lon=as.numeric(gsub(lonlat_pat, "\\1", x)), lat=as.numeric(gsub(lonlat_pat, "\\2", x)))
	# from_to[,lon:=as.numeric(gsub(lonlat_pat, "\\1", from_site))]
	# from_to[,lat:=as.numeric(gsub(lonlat_pat, "\\2", from_site))]
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

siteLL <- from_to[,stratum2ll(unique(c(from_site,to_site))),by='reg']
fromLL <- from_to[,stratum2ll(from_site)[,list(from_lon=lon,from_lat=lat)]]
toLL <- from_to[,stratum2ll(to_site)[,list(to_lon=lon,to_lat=lat)]]
from_to <- cbind(from_to, fromLL, toLL)
from_to[,dist:=geoDist(x=data.table(from_lon, from_lat), y=data.table(to_lon, to_lat))]
# from_to[,.SD[which.min(dist)],by=c('reg','spp','from_year','to_year')] # if multiple have smallest, only returns first
# from_to <- from_to[,.SD[dist%in%min(dist)],by=c('reg','spp','from_year','to_year')] # if multiple have smallest, returns all of them

getNB <- function(lon, lat){
	locs <- simplify2array(ll2km(lon, lat))[,2:1] # converts to km while accounting for lon-lat coordinates
	nn1 <- spdep::knn2nb(spdep::knearneigh(locs, k=1))
	max2NDist <- max(unlist(spdep::nbdists(nn1, locs)))
	localAC$max2NDist <- max2NDist
	localAC$nb <- spdep::dnearneigh(locs, d1=0, d2=max2NDist) # graph2nb(gabrielneigh(locs))
}

eval(figure_setup())
ur <- from_to[,names(pretty_reg)[names(pretty_reg)%in%unique(reg)]]
nbList <- structure(vector("list", length(ur)), .Names=ur)
for(r in 1:length(ur)){
	nbList[[ur[r]]] <- siteLL[reg==ur[r],getNB(lon, lat)]
}


eval(figure_setup())
ur <- from_to[,names(pretty_reg)[names(pretty_reg)%in%unique(reg)]]
g_list <- structure(vector("list", 9), .Names=ur)
for(r in 1:length(ur)){
	tr <- ur[r]
	tdat <- from_to[reg==tr,list(from=from_site, to=to_site, reg, spp, from_year, to_year)]
	tsll <- siteLL[reg==tr][,reg:=NULL]
	vertices <- tsll[tdat[,unique(c(from,to))], on="site"]
	tg <- igraph::graph_from_data_frame(tdat, directed=TRUE, vertices=vertices)
	g_list[[tr]]$graph <- tg
	g_list[[tr]]$layout <- as.matrix(vertices[,list(x=lon,y=lat)])
}
# plot(mapOwin[['neus']])
# map(add=TRUE, fill=TRUE, col='gray')
# qgraph(as_adj(g_list[['neus']]$graph), layout=g_list[['neus']]$layout, labels=FALSE, rescale=FALSE, plot=FALSE, color=NA, arrows=3, threshold=2, border.color=adjustcolor('gray',0.25))
# write.csv(as.matrix(as_adj(g_list[['neus']]$graph)), file="../manuscript/adjMat_example_NortheastUS.csv", row.names=TRUE)
# write.csv(as.matrix(as_adj(g_list[['ebs']]$graph)), file="../manuscript/adjMat_example_EasternBeringSea.csv", row.names=TRUE)
# write.csv(as.matrix(as_adj(g_list[['ai']]$graph)), file="../manuscript/adjMat_example_AleutianIslands.csv", row.names=TRUE)

dev.new(width=7, height=3)
eval(figure_setup())
# ur <- from_to[,unique(reg)[unique(reg)%in%names(pretty_reg)]]
ur <- from_to[,names(pretty_reg)[names(pretty_reg)%in%unique(reg)]]
map_layout <- trawl_layout()
par(mar=c(0.9,0.9,0.25,0.25), mgp=c(0.5,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.95,0.7,1,0.1))
layout(map_layout)
# par(mfrow=c(3,3), mar=c(0.5,0.5,0.5,0.5))
for(r in 1:length(ur)){
	tr <- ur[r]
	plot(mapOwin[[tr]], main="")
	maps::map(add=TRUE, fill=TRUE, col='lightgray')
	am <- as_adj(g_list[[tr]]$graph)
	
	ecol <- as.matrix(am)
	ecol[] <- 'blue'
	sll <- gsub(" [0-9]{1,4}$", "", rownames(ecol))
	isSame <- outer(sll, sll, "==")
	ecol[isSame] <- 'red'
	
	thresh <- trunc(median(unique(as.numeric(am))))
	qgraph(am, layout=g_list[[tr]]$layout, labels=FALSE, rescale=FALSE, plot=FALSE, normalize=TRUE, vTrans=0, diag=TRUE, edge.width=1, loop=30, asize=5, directed=TRUE, arrows=2, edge.color=ecol, border.color=adjustcolor('gray',0.25), threshold=0, vsize=10)
}


# ==========
# = Banner =
# ==========


eval(figure_setup())
# ur <- from_to[,unique(reg)[unique(reg)%in%names(pretty_reg)]]
ur <- from_to[,names(pretty_reg)[names(pretty_reg)%in%unique(reg)]]
map_layout <- trawl_layout()
par(mar=c(0.9,0.9,0.25,0.25), mgp=c(0.5,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.95,0.7,1,0.1))
layout(map_layout)
# par(mfrow=c(3,3), mar=c(0.5,0.5,0.5,0.5))
for(r in 1:length(ur)){
	tr <- ur[r]
	plot(mapOwin[[tr]], main="")
	maps::map(add=TRUE, fill=TRUE, col='lightgray')
	mtext(paste(dim(as_adj(g_list[[tr]]$graph)), collapse=" by "), side=3, line=-1, font=2, cex=2)
}




# =================
# = patrick magic =
# =================
rowScale <- function(x){
	sf <- function(y){
		y/pmax(sum(y),1)
	}
	o <- t(apply(x, 1, sf))
	stopifnot(all(rowSums(o)==1 | rowSums(o)==0))
	return(o)
}
neus_pat <- as.matrix(as_adj(g_list[['neus']]$graph))
neus_patScale <- rowScale(neus_pat) #t(t(neus_pat)/pmax(rowSums(neus_pat),1))
which(neus_patScale%^%1E5 != 0)
all(eigen(neus_patScale)[[2]][,1] == 0)

