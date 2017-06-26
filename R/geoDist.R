#' Geographic Distance
#' Find the geographic distance between two sets of coordinates
#' @param x lon
#' @param y lat
#' 
#' @return a numeric
#' @export
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


#' Convert Stratum to Lon-Lat
#' Convert a stratum name (in my format) to longitude and latitude by regex
#' 
#' @param x a character vector of stratum names
#' 
#' @return a data.table with columns for site (stratum), lon, lat
#' 
#' @export
stratum2ll <- function(x){
	lonlat_pat <- "^(-[0-9]{2,3}\\.[0-9]{2}) ([0-9]{2,3}\\.[0-9]{2}) [0-9]{1,4}$"
	data.table(site=x, lon=as.numeric(gsub(lonlat_pat, "\\1", x)), lat=as.numeric(gsub(lonlat_pat, "\\2", x)))
	# from_to[,lon:=as.numeric(gsub(lonlat_pat, "\\1", from_site))]
	# from_to[,lat:=as.numeric(gsub(lonlat_pat, "\\2", from_site))]
}