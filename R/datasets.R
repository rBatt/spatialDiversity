

#' Data from All Regions
#' 
#' Data from all regions generated by \code{get_data_all}
#' 
#' @format data.table
"data_all2"

#' Data for Maps
#' 
#' Spatial information related to colonization, extinction, richness, depth, and bottom temperature; does not contain any explicit temporal information
#' 
#' @details
#' Meant for making maps. Other processed data sets emphasize the temporal dimension. This data set collapses time, and emphasizes 2D spatial information.
#' 
#' \tabular{lll}{
#' \code{stratum} \tab character \tab sampling stratum, redefined via \code{\link{trim_msom}} in \code{\link{get_data_all}} (and therefore in \code{\link{data_all}}) \cr
#' \code{lon} \tab numeric \tab longitude in degrees East (such that negative values are in the Western hemisphere) \cr
#' \code{lat} \tab numeric \tab latitude \cr
#' \code{depth} \tab numeric \tab depth in meters; in each year hauls in the stratum are averaged, then this value is the cross-year average of those annual depths \cr
#' \code{n_spp_col_weighted, n_spp_ext_weighted} \tab numeric \tab total colonizers per stratum colonized (# species colonizing this stratum / # of strata colonized by each species) in this stratum, averaged across years; similar for ext= extinction \cr
#' \code{reg} \tab character \tab region \cr
#' \code{avgRich, sdRich} \tab numeric \tab cross-year average (standard deviation) richness in this stratum, based on naive richness (not MSOM estimate) \cr
#' \code{avgBtemp, sdBtemp} \tab numeric \tab cross-year average (standard deviation) of bottom temperature in this stratum \cr
#' }
#' 
#' @seealso \code{\link{spp_master}}, \code{\link{comm_master}}
#' 
#' @format data.table
"mapDat"

#' Regional Polygon Outlines
#' 
#' Like a convex hull, but a bit more form-fitting. They should all be counterclockwise (to work with \code{spatstat::ppp}), and they are not quite closed. These were manually created.
#' 
#' @details
#' Meant for making maps. Other processed data sets emphasize the temporal dimension. This data set collapses time, and emphasizes 2D spatial information.
#' 
#' @details
#' \tabular{lll}{
#' \code{reg} \tab character \tab region \cr
#' \code{lonP} \tab numeric \tab longitude in degrees East (such that negative values are in the Western hemisphere) \cr
#' \code{latP} \tab numeric \tab latitude \cr
#' }
#' 
#' @format data.table
"outlines"

#' Region outlines in owin format
#' 
#' Transformed \code{\link{outlines}} for use with the \code{spatstat package}
#' 
#' @seealso \code{\link{make_owin}}
#' 
#' @format owin
"mapOwin"

#' Results for local autocorrelation
#' 
#' Includes distance, neighborhood, and local autocorrelation (and associated pvalues, z scores, lon lat) for each site
#' 
#' @seealso \code{\link{spatial_ac}}
#' 
#' @format list data.table
"localAC"

#' Timing of species' colonizations and extinctions
#' 
#' Includes species, year, region, and whether it colonized or went extinct (or both)
#' 
#' @seealso \code{\link{data_all2}} \code{\link{mapDat}}
#' 
#' @format list data.table
"col_ext_dt"

#' A lightweight equivalent of trawlDiversity's spp_master
#' 
#' Includes species, year, region, and the colonization/ extinction events
#' 
#' @format list data.table
"spp_master2"

#' Colonization and Extinction Sites, Distances
#' 
#' The sites (strata) of colonizations, extinctions, and the distances between these pairs of sites. Colonization-extinction pairs are on a per-species basis (and per region), extinction events always follow colonization events, the species has to be present for more than 1 year in a row. Note that only species that both colonize and go extinct are included here.
#' 
#' @format list data.table
"from_to"





