
# ---- Map Data ----
mapDat <- make_mapDat(p)

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
mapOwin <- trawlDiversity::make_owin(mapDat, outlines)

# ---- calculate spatial autocorrelation ----
rs <- mapDat[,una(reg)]
nr <- length(rs)
localAC <- list(colonization=list(), extinction=list())
# lac_val <- c("n_spp_col_weighted", "n_spp_col_unique")[1]
lac_val <- c("n_spp_col_weighted", "n_spp_ext_weighted")
ce_types <- names(localAC) #c("colonization","extinction")
for(ce in 1:2){
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


