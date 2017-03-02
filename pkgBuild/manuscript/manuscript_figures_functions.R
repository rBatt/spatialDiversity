
library('maps')
library('raster')
library('spatstat')
library("spdep")
library('rbLib')
library('spatialDiversity')
library("data.table")

cePCH <- function(region, ce_type=c("tot","u","original")){
	ce_type <- match.arg(ce_type)
	if(ce_type=="original"){
		ce_type <- ""
		ce_suff <- c("col","ext")
	}else{
		ce_suff <- paste0(ce_type, c("Col","Ext"))
	}
	
	mD <- mapDat[reg==region]
	setnames(mD,c("n_spp_col_weighted", "n_spp_ext_weighted"), c("col", "ext"))
	
	pval_col <- mD[,get(paste0("lI_pvalue_",ce_suff[1]))]
	pval_ext <- mD[,get(paste0("lI_pvalue_",ce_suff[2]))]
	sigColInd <- pval_col<0.05
	sigExtInd <- pval_ext<0.05
	
	muCol <- mD[,mean(get(ce_suff[1]))]
	muExt <- mD[,mean(get(ce_suff[2]))]
	
	hotspotIndCol <- sigColInd #& (totCol > muCol)
	hotspotIndExt <- sigExtInd #& (totExt > muExt)
	both <- hotspotIndExt&hotspotIndCol
	neither <- !hotspotIndExt&!hotspotIndCol
	
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
	
	return(list(pchs, cols))
}



#
#
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

# ---- Colonization Rate ----
ceRate_map <- function(ce=c("colonization","extinction","richness","uCol","uExt","totCol","totExt"), add_lisa=TRUE, ce_type=c("tot","u","original")){
	ce <- match.arg(ce)
	eval(figure_setup())
	map_layout <- trawl_layout()
	par(mar=c(0.9,0.9,0.25,0.25), mgp=c(0.5,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.85,0.6,1,0.1))
	layout(map_layout)
	
	ce_type <- match.arg(ce_type) 
	ce_type <- switch(ce,
		colonization = 'original',
		extinction = 'original',
		richness = ce_type,
		uCol = 'u',
		uExt = 'u',
		totCol = 'tot',
		totExt = 'tot'
	)

	map_col <- grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	toRast <- function(p){
		crs_orig <- CRS("+proj=longlat")
		
		x0 = rep(p$xcol, length(p$yrow))
		y0 = rep(p$yrow, each=length(p$xcol))
		# xy <- mapproject(x0, y0, proj="albers", parameters=extent(r)[3:4], orientation=c(mean(y0), mean(x0), 0))
		z = c(t(as.matrix(p)))
		
		r0 <- raster::rasterFromXYZ(cbind(x0, y0, z), crs=crs_orig)
		
		# crs_new <- CRS(paste0("+proj=aea +lat_1=", min(y0), " +lat_2=", max(y0), " +lat_0=", mean(y0), " +lon_0=", mean(x0), " +x_0=0 +y_0=0 +ellps=WGS84 +datum=NAD83 +units=m +no_defs"))
		# proj_to <- projectExtent(r0, crs=crs_new)
		# r <- projectRaster(from=r0, to=proj_to)
		# plot(r, asp=1)
		# map(add=T, proj="albers", parameters=extent(r)[3:4], orientation=c(mean(y0), mean(x0), 0))
		r <- r0
		return(r)
	}

	u_regs <- mapDat[,unique(reg)]
	rs <- mapDat[,una(reg)]
	nr <- length(rs)
	mapPPP_ce <- list()
	for(r in 1:nr){
		tr <- rs[r]
		td <- mapDat[reg==tr]
		if(ce=="colonization"){
			mapPPP_ce[[tr]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_col_weighted], window=mapOwin[[tr]]) 
		}else if(ce=="extinction"){
			mapPPP_ce[[tr]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_ext_weighted], window=mapOwin[[tr]])
		}else if(ce=="richness"){
			mapPPP_ce[[tr]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,avgRich], window=mapOwin[[tr]])
		}else{
			mapPPP_ce[[tr]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,get(ce)], window=mapOwin[[tr]])
		}
		
		t_idw <- spatstat::Smooth(mapPPP_ce[[tr]], hmax=1)
		z <- toRast(t_idw)
		raster::image(z, col=map_col, xlab="", ylab="", asp=1)
		if(tr=="ai"){
			w2 <- map('world2', plot=FALSE)
			w2$x <- w2$x-360
			map(w2, add=TRUE, fill=TRUE, col="lightgray")
		}else{
			map(add=TRUE, fill=TRUE, col="lightgray")
		}
		
		if(add_lisa){
			pc <- cePCH(tr, ce_type=ce_type)
			pc[[1]][pc[[1]]==1] <- NA
			
			sigRichInd <- mapDat[reg==tr, lI_pvalue_rich<0.05]
			hotspotIndRich <- sigRichInd
			
			mapDat[reg==tr, points(lon, lat, pch=pc[[1]], col="white", lwd=3, cex=1.2)]
			mapDat[reg==tr,points(lon[hotspotIndRich], lat[hotspotIndRich], col='white', pch=20, cex=1.5)]
			
			mapDat[reg==tr, points(lon, lat, pch=pc[[1]], col=pc[[2]], cex=1.2, lwd=1)]
			mapDat[reg==tr,points(lon[hotspotIndRich], lat[hotspotIndRich], col='white', pch=20, cex=1)]
			mapDat[reg==tr,points(lon[hotspotIndRich], lat[hotspotIndRich], col='black', pch=20, cex=0.7)]
		}
		
	
		zl <- range(values(z)*10, na.rm=TRUE)
		switch(tr,
			ebs = mapLegend(x=0.025, y=0.225, h=0.37, w=0.025, zlim=zl, cols=map_col, lab.cex=1),
			ai = mapLegend(x=0.985, y=0.3, w=0.02, h=0.55, zlim=zl, cols=map_col, lab.cex=1),
			goa = mapLegend(x=0.985, y=0.15, w=0.02,  zlim=zl, cols=map_col, lab.cex=1),
			wctri = mapLegend(x=0.1, y=0.1125, w=0.10, h=0.175, zlim=zl, cols=map_col, lab.cex=1),
			gmex = mapLegend(x=0.95, y=0.225, h=0.375, zlim=zl, cols=map_col, lab.cex=1),
			sa = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			neus = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			shelf = mapLegend(x=0.95, y=0.185, zlim=zl, h=0.31, cols=map_col, lab.cex=1),
			newf = mapLegend(x=0.05, y=0.125, h=0.20, zlim=zl, cols=map_col, lab.cex=1)
		)
		switch(tr,
			ebs = legend("topright", legend="A", bty='n', text.font=2, inset=c(-0.02,-0.15), cex=1.25, text.col='black'),
			ai = legend("topleft", legend="C", bty='n', text.font=2, inset=c(-0.065,-0.2), cex=1.25, xpd=T),
			goa = legend("topleft", legend="B", bty='n', text.font=2, inset=c(-0.065,-0.06), cex=1.25),
			wctri = legend("topright", legend="E", bty='n', text.font=2, inset=c(-0.01,-0.05), cex=1.25, text.col='black'),
			gmex = legend("topleft", legend="G", bty='n', text.font=2, inset=c(-0.175,-0.12), cex=1.25, text.col='black'),
			sa = legend("topleft", legend="I", bty='n', text.font=2, inset=c(-0.15,-0.075), cex=1.25, text.col='black'),
			neus = legend("topleft", legend="H", bty='n', text.font=2, inset=c(-0.125,-0.05), cex=1.25, text.col='black'),
			shelf = legend("topleft", legend="D", bty='n', text.font=2, inset=c(-0.1,-0.125), cex=1.25, text.col='black'),
			newf = legend("topright", legend="F", bty='n', text.font=2, inset=c(-0.01,-0.05), cex=1.25)
		)
	}
	if(ce=="colonization"){
		mtext(bquote(Colonization~Rate~(C[w]~~decade^-1)), side=3, outer=TRUE, font=2, line=-0.3)
	}else if(ce=="extinction"){
		mtext(bquote(Extinction~Rate~(E[w]~~decade^-1)), side=3, outer=TRUE, font=2, line=-0.3)
	}else if(ce=="richness"){
		mtext(bquote(Observed~Richness), side=3, outer=TRUE, font=2, line=-0.3)
	}else{
		mtext(ce, side=3, outer=TRUE, font=2, line=-0.3)
	}
	mtext(bquote(Longitude~(phantom()*degree*E)), side=1, line=0.15, outer=TRUE)
	mtext(bquote(Latitude~(phantom()*degree*N)), side=2, line=-0.4, outer=TRUE)
	invisible(NULL)
}

# ---- Neighborhood used in Moran's I ----
nb_moranI <- function(ce=c("richness", "colonization", "extinction", "uCol", "uExt", "totCol", "totExt")){
	eval(figure_setup())
	map_layout <- trawl_layout()
	par(mar=c(0.25,0.25,0.25,0.25), mgp=c(0.25,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.1,0.1,0.1,0.1))
	layout(map_layout)
	pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West\nCoast\nUS", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")

	rs <- mapDat[,una(reg)]
	nr <- length(rs)
	for(r in 1:nr){
		tr <- rs[r]
		t_lac <- localAC[[ce]][[tr]]
		plot(mapOwin[[tr]], coords=t_lac$I[,list(lon,lat)], add=FALSE, main="")
		box()
		if(rs[r]=='wctri'){
			mtext(pretty_reg[tr], side=3, line=-3)
		}else{
			mtext(pretty_reg[tr], side=3, line=-1)
		}
		plot(t_lac$nb, t_lac$I[,list(lon,lat)], add=TRUE, col=adjustcolor('black', 0.5), cex=0.8, lwd=0.5)
		sig_lac <- t_lac$I[,lI_pvalue]<0.05
		if(any(sig_lac)){
			zl <- range(t_lac$I[sig_lac,Ii], na.rm=TRUE)
			t_col <- rbLib::zCol(256,t_lac$I[sig_lac,Ii])
			map_col <- rbLib::zCol(6, 1:6)
			points(x=t_lac$I[sig_lac,lon], y=t_lac$I[sig_lac,lat], bg=t_col, pch=21, cex=1.1)
			switch(tr,
				ebs = mapLegend(x=0.025, y=0.225, h=0.37, w=0.025, zlim=zl, cols=map_col, lab.cex=1),
				ai = mapLegend(x=0.985, y=0.3, w=0.02, h=0.55, zlim=zl, cols=map_col, lab.cex=1),
				goa = mapLegend(x=0.985, y=0.15, w=0.02,  zlim=zl, cols=map_col, lab.cex=1),
				wctri = mapLegend(x=0.1, y=0.1125, w=0.10, h=0.175, zlim=zl, cols=map_col, lab.cex=1),
				gmex = mapLegend(x=0.95, y=0.225, h=0.375, zlim=zl, cols=map_col, lab.cex=1),
				sa = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
				neus = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
				shelf = mapLegend(x=0.95, y=0.185, zlim=zl, h=0.31, cols=map_col, lab.cex=1),
				newf = mapLegend(x=0.05, y=0.125, h=0.20, zlim=zl, cols=map_col, lab.cex=1)
			)
		}else{
			# plot(x=locs[,1], y=locs[,2], col='blue')
		}
	}
	invisible(NULL)
}


