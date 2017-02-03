
update_dependencies <- function(){
	devtools::use_package("data.table", type="Depends") # Basis for handling all data sets
	devtools::use_package("trawlData", type="Depends") # Meets basic requirements for data content and format
	devtools::use_package("trawlDiversity", type="Depends") # builds off of previous analysis; uses some data objects
	
	devtools::use_package("rbLib", type="Imports")
	devtools::use_package("stats", type="Imports")
	devtools::use_package("methods", type="Imports")	
}