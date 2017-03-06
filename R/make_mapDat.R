#' Make data for maps related to diversity analysis
#' 
#' For each region, get the data appropriate for maps of richness and environmental metrics.
#' 
#' @param p The \code{p} data object from the process_msom_script.R file.
#' 
#' @details
#' The data.table contains summary statistics to collapse time (long-term mean, e.g.), and retains stratum-specific information.
#' 
#' @return a data.table
#' 
#' @export
make_mapDat <- function(p){
	mapDat <- list()
	for(r in 1:length(p)){
		# part 1
		pt1 <- trawlAgg(
			p[[r]]$rd,
			bioFun=meanna, envFun=meanna,
			envCols=c("btemp","depth","stemp","lon","lat"),
			bio_lvl="spp",time_lvl="year",space_lvl="stratum",
			metaCols=c("reg"),meta.action="unique1"	
		)
		
		# part 2
		pt2 <- pt1[,j={
			avgRich <- .SD[,trawlData::lu(spp),by="time_lvl"][,meanna(V1)]
			sdRich <- .SD[,trawlData::lu(spp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			avgBtemp <- .SD[,meanna(btemp),by="time_lvl"][,meanna(V1)]
			sdBtemp <- .SD[,meanna(btemp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			data.table(avgRich=avgRich, sdRich=sdRich, avgBtemp=avgBtemp, sdBtemp=sdBtemp)
		},by=c("reg","stratum")]
		
		# part 3, unique colonizers (belongs in trawlDiversity::get_colonizers, if put in this pkg)
		n_uSppPerStrat_expr <- bquote(.SD[(col_logic),length(unique(spp))])
		get_uniqueCE <- function(x){
			reg_name <- x$rd[,unique(reg)]
			t_c <- x$colonization
			n_uSppPerStrat_exts <- t_c$ext_dt[, list(uExt=as.numeric(eval(n_uSppPerStrat_expr))), by=c("stratum")]
			n_uSppPerStrat_cols <- t_c$col_dt[, list(uCol=as.numeric(eval(n_uSppPerStrat_expr))), by=c("stratum")]
			data.table(reg=reg_name, merge(n_uSppPerStrat_exts,n_uSppPerStrat_cols, all=TRUE))
		}
		n_uSppPerStrat_ce <- get_uniqueCE(p[[r]]) # n_uSppPerStrat_ce <- rbindlist(lapply(p, get_uniqueCE))
		
		# part 4, total colonizers (belongs in trawlDiversity::get_colonizers, if put in this pkg)
		n_totSppPerStrat_expr <- bquote(.SD[,sum(col_logic)])
		get_totalCE <- function(x){
			reg_name <- x$rd[,unique(reg)]
			t_c <- x$colonization
			n_totSppPerStrat_exts <- t_c$ext_dt[, list(totExt=as.numeric(eval(n_totSppPerStrat_expr))), by=c("stratum")]
			n_totSppPerStrat_cols <- t_c$col_dt[, list(totCol=as.numeric(eval(n_totSppPerStrat_expr))), by=c("stratum")]
			data.table(reg=reg_name, merge(n_totSppPerStrat_exts,n_totSppPerStrat_cols, all=TRUE))
		}
		n_totSppPerStrat_ce <- get_totalCE(p[[r]])
		
		
		# Merge parts
		to_merge <- c(p[[r]]$colonization[c("n_spp_col_weighted_tot","n_spp_ext_weighted_tot","n_cep")])
		mapDat[[r]] <- merge(to_merge[["n_spp_col_weighted_tot"]], to_merge[["n_spp_ext_weighted_tot"]], by=c("stratum","lon","lat","depth","yrs_sampled"),all=TRUE)
		mapDat[[r]] <- merge(mapDat[[r]], pt2, by="stratum",all=TRUE)
		mapDat[[r]] <- merge(mapDat[[r]], n_uSppPerStrat_ce, by=c("reg","stratum"),all=TRUE)
		mapDat[[r]] <- merge(mapDat[[r]], n_totSppPerStrat_ce, by=c("reg","stratum"),all=TRUE)
	}
	mapDat <- rbindlist(mapDat)[reg!="wcann"]
	# mapDim <- mapDat[,list(r_lon=diff(range(lon)),r_lat=diff(range(lat))),by="reg"]
	# mapDim[,c("lon_scale","lat_scale"):=list(r_lon/min(r_lon), r_lat/min(r_lat))]
	# mapDim[,ll_ratio:=r_lon/r_lat]
	
	# put u/tot values in units of "per year" to control for some sites not being sampled in all years
	mapDat[,uCol:=uCol/yrs_sampled,by="reg"]
	mapDat[,uExt:=uExt/yrs_sampled,by="reg"]
	mapDat[,totCol:=totCol/yrs_sampled,by="reg"]
	mapDat[,totExt:=totExt/yrs_sampled,by="reg"]
	
	return(mapDat)
}