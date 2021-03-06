#' Make data for annual maps 
#' 
#' Similar as for mapDat, but with a column for time
#' 
#' @param p The \code{p} data object from the process_msom_script.R file.
#' 
#' @details
#' The data.table contains spatial and temporal information on colonizations, extinctions, richness, and the species involved. Some environmental information, too.
#' 
#' @return a data.table
#' 
#' @export
make_mapDatTime <- function(p){
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
			list(rich=trawlData::lu(spp))
		},by=c("reg","stratum","time_lvl")]
		
		# part 3, look at the strata for colonizations and extinctions in each year
		annCE_expr <- bquote(.SD[(col_logic),length(unique(spp))])
		get_annCE <- function(x){
			reg_name <- x$rd[,unique(reg)]
			t_c <- x$colonization
			annExts <- t_c$ext_dt[, list(uExt=eval(annCE_expr)), by=c("stratum","time_lvl")]
			annCols <- t_c$col_dt[, list(uCol=eval(annCE_expr)), by=c("stratum","time_lvl")]
			data.table(reg=reg_name, merge(annExts,annCols, all=TRUE))
		}
		annCE <- get_annCE(p[[r]]) # n_uSppPerStrat_ce <- rbindlist(lapply(p, get_uniqueCE))
		
		
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
	
	return(mapDat)
}