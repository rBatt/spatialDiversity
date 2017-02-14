
library(trawlData)
library(trawlDiversity)
library(spatialDiversity)
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

mapDat <- spatialDiversity::make_mapDat(small_p)

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

# =================================
# = Save Data Ojbects for Package =
# =================================
save(mapDat, file="data/mapDat.RData")
save(mapOwin, file="data/mapOwin.RData")
save(localAC, file="data/localAC.RData")


