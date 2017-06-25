
library(trawlData)
library(trawlDiversity)
# library(spatialDiversity)
library(rbLib)

setwd("~/Documents/School&Work/pinskyPost/spatialDiversity")
invisible(lapply(list.files("R",full=TRUE), source))

# ---- Map Data ----
u_dat <- trawlDiversity::comm_master[,list(msom_years=unique(year)), by="reg"]
u_reg <- u_dat[,unique(reg)]
n_reg <- length(u_reg)
small_p <- vector("list", n_reg)
data_all2 <- copy(trawlDiversity::data_all)
data_all2 <- data_all2[K==1]

# par(mfrow=c(3,3))
# data_all2[,unique(Kmax),by=c('reg','year','stratum')][,sum(V1<=1)/length(V1),by=c("reg","year")][reg!='wcann',j={plot(year,V1,ylim=c(0,1));NULL},by='reg']
# data_all2[,unique(Kmax),by=c('reg','year','stratum')][,sum(V1<=1)/length(V1),by=c("reg","year")][reg!='wcann',mean(V1),by='reg']
# # from this i learn that in *almost* every year and region there is at least 1 site with only 1 tow. There are a couple instances where the minimum number of tows might be a bit higher (like 2). The fraction of sites with 1 tow is around 20%; averaging across years, some regions typicall have as few as 7% of sites with <= 1 tow (gmex, goa at 8%), or with a cross-year average of as many as 47% of sites having <= 1 tow (ebs; next in line is newf with 29%, so ebs is an extreme case).

for(r in 1:n_reg){
	t_reg <- u_reg[r]
	small_p[[r]] <- trawlDiversity::process_obsRich(data_all2[reg==t_reg], msom_yrs=u_dat[reg==t_reg,msom_years])
}

col_ext_dt <- lapply(small_p, function(x)data.table(reg=x$rd[,unique(reg)], x$colonization$col_ext_dt))
col_ext_dt <- rbindlist(col_ext_dt)

# # start new
# n_uSppPerStrat_expr <- bquote(.SD[(col_logic),length(unique(spp))])
# get_uniqueCE <- function(x){
# 	reg_name <- x$rd[,unique(reg)]
# 	t_c <- x$colonization
# 	n_uSppPerStrat_exts <- t_c$ext_dt[, list(uExt=eval(n_uSppPerStrat_expr)), by=c("stratum")]
# 	n_uSppPerStrat_cols <- t_c$col_dt[, list(uCol=eval(n_uSppPerStrat_expr)), by=c("stratum")]
# 	data.table(reg=reg_name, merge(n_uSppPerStrat_exts,n_uSppPerStrat_cols, all=TRUE))
# }
# n_uSppPerStrat_ce <- rbindlist(lapply(small_p, get_uniqueCE))
#
# # # quick plot
# # n_uSppPerStrat_ce[,plot(jitter(uExt, factor=2), jitter(uCol, factor=2), col=as.factor(reg), pch=20)]
# # abline(a=0, b=1)
# # n_uSppPerStrat_ce[,legend("topleft", legend=unique(reg), col=as.factor(unique(reg)), pch=19)]
# #
# # par(mfrow=c(3,3))
# # n_uSppPerStrat_ce[,j={plot(uExt, uCol);abline(a=0,b=1);mtext(reg, side=3, line=0.5, font=2)},by='reg']
#
# # end new

mapDat <- make_mapDat(small_p) # spatialDiversity::

# # start new
# mapDat2 <- merge(mapDat, n_uSppPerStrat_ce)
#
# par(mfrow=c(3,3))
# mapDat2[,j={plot(uCol,n_spp_col_weighted*yrs_sampled);abline(a=0,b=1);mtext(reg, side=3, line=0.5, font=2)},by="reg"]
#
# par(mfrow=c(3,3))
# mapDat2[,j={plot(uExt,n_spp_ext_weighted*yrs_sampled);abline(a=0,b=1);mtext(reg, side=3, line=0.5, font=2)},by="reg"]
# # end new




# ---- Function to help outline a region ----
regOutline <- function(X){
	dev.new()
	outlines <- list()
	rs <- X[,una(reg)]
	nr <- length(rs)
	for(r in 1:nr){
		td <- X[reg==rs[r]]
		rast <- raster(xmn=td[,min(lon)-1], xmx=td[,max(lon)+1], ymn=td[,min(lat)-1], ymx=td[,max(lat)+1], res=c(0.5,0.5))
		td[,plot(rasterize(cbind(x=lon,y=lat), rast))]
		td[,points(lon, lat, cex=5)]
		map(add=TRUE, fill=F)
		to <- locator(type='o', pch=20, col='blue')
		outlines[[r]] <- data.table(reg=rs[r], lonP=to$x, latP=to$y)
	}
	return(rbindlist(outlines))
}
# outlines <- regOutline(mapDat)
# save(outlines, file="~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/data/outlines.RData")

# ---- make window polygon for ppp ----
mapOwin <- make_owin(mapDat, outlines)

# ---- calculate spatial autocorrelation ----
rs <- mapDat[,una(reg)]
nr <- length(rs)
localAC <- list(richness=list(), colonization=list(), extinction=list(), uCol=list(), uExt=list(), totCol=list(), totExt=list())
# lac_val <- c("n_spp_col_weighted", "n_spp_col_unique")[1]
lac_val <- c("avgRich", "n_spp_col_weighted", "n_spp_ext_weighted", "uCol", "uExt", "totCol", "totExt")
ce_types <- names(localAC) #c("colonization","extinction")
for(ce in 1:length(lac_val)){
	for(r in 1:nr){
		lac_val_ce <- lac_val[ce]
		# t_lac <- with(mapDat[reg==rs[r]], spatial_ac(lon, lat, n_spp_col_weighted))
		t_lac <- with(mapDat[reg==rs[r]][complete.cases(mapDat[reg==rs[r]])], spatial_ac(lon, lat, eval(s2c(lac_val_ce))[[1]]))
		t_lac$I <- data.table(mapDat[reg==rs[r], list(reg,stratum)], t_lac$I)
		localAC[[ce_types[ce]]][[rs[r]]] <- t_lac
	}
	lac_2mapDat <- rbindlist(lapply(localAC[[ce_types[ce]]], function(x1)x1$I))
	
	if(ce_types[ce]=="colonization"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_col=Ii,lI_pvalue_col=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="extinction"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_ext=Ii,lI_pvalue_ext=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="richness"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_rich=Ii,lI_pvalue_rich=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="uCol"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_uCol=Ii,lI_pvalue_uCol=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="uExt"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_uExt=Ii,lI_pvalue_uExt=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="totCol"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_totCol=Ii,lI_pvalue_totCol=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}else if(ce_types[ce]=="totExt"){
		mapDat <- merge(mapDat, lac_2mapDat[,list(reg,stratum,Ii_totExt=Ii,lI_pvalue_totExt=lI_pvalue)], by=c("reg","stratum"), all=TRUE)
	}
	mapDat[,reg:=factor(reg, levels=c("ebs", "ai", "goa", "wctri", "gmex", "sa", "neus", "shelf", "newf"))]
	setorder(mapDat, reg, stratum)
	mapDat[,reg:=as.character(reg)]
}

# ===============
# = spp_master2 =
# ===============
spp_master2 <- merge(data_all2, col_ext_dt, by=c("reg","year","spp"), all=TRUE)
spp_master2[is.na(col) & is.na(ext), c('col','ext'):=list(0, 0)] # for the 'neither' spp
define_ce_categ <- function(X){
	ext <- X[,ext]
	col <- X[,col]
	if(all(col==0) & all(ext==0)){
		return("neither")
	}
	if(all(col==0) & any(ext==1)){
		return("leaver")
	}
	if(any(col==1) & all(ext==0)){
		return("colonizer")
	}
	if(any(col==1) & any(ext==1)){
		return("both")
	}
}
spp_master2[,ce_categ:=define_ce_categ(.SD),by=c("reg","spp")]


# =============================================================
# = Sites of Colonization (from) and Extinction (to): from_to =
# =============================================================
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

from_to_years <- spp_master2[ce_categ=="both",j={
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
pt1 <- unique(spatialDiversity::data_all2[,list(reg,spp,year,stratum)])[from_to_years,on=c("reg","spp",year="from_year")]
pt1 <- pt1[,list(reg,spp,from_year=year,to_year,from_site=stratum)]
pt2 <- unique(spatialDiversity::data_all2[,list(reg,spp,year,stratum)])[from_to_years,on=c("reg","spp",year="to_year")]
pt2 <- pt2[,list(reg,spp,from_year,to_year=year,to_site=stratum)]
from_to <- merge(pt1, pt2, all=TRUE, allow.cartesian=TRUE)
from_to <- from_to[!from_year==to_year] # remove when col ext happen same year

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


# =================================
# = Save Data Ojbects for Package =
# =================================
save(mapDat, file="data/mapDat.RData")
save(mapOwin, file="data/mapOwin.RData")
save(localAC, file="data/localAC.RData")
save(data_all2, file="data/data_all2.RData")
save(col_ext_dt, file="data/col_ext_dt.RData")
save(spp_master2, file="data/spp_master2.RData")
save(from_to, file="data/from_to.RData")


