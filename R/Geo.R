### COLLECTION OF GEOGRAPHIC INFORMATION FUNCTIONS ###

#.coords2country =====
#' @title Extracts country ISO3 based on locations
#' @param xydat \emph{data.frame} with x and y coordinates
#' @param .countries.shapefile \emph{multipolygon} sf objects  of world countries and their associated ISO3 codes
#' @param .points.crs \emph{character} crs for the occurrence data
#' @param xf \emph{character} name of the x coordinate field
#' @param yf \emph{character} name of the y coordinate field
#' @param ctryNameField \emph{character} Column name in .countries shapefile where ISO3 are indicated
#' @param verbose \emph{logical} Print messages?
#' @return \emph{character} factor vector with ISO3 codes for countries 
#' @family Geo
#' @author Josep M Serra Diaz


.coords2country = function(xydat,
                           xf=xf,
                           yf=yf,
                           .countries.shapefile=NULL,
                           .points.crs=NULL,
                           ctryNameField=NULL,verbose=FALSE){  
  if (is.null(.countries.shapefile)) {   .countries.shapefile = sf::st_as_sf(rnaturalearthdata::countries110) ;  ctryNameField ='ISO_A3'}
  if (!inherits(.countries.shapefile,'sf')) {stop (".countries shapefile not loaded ")}
  if (is.null(.points.crs)) {.points.crs <- sf::st_crs(.countries.shapefile); if(verbose) print (paste ('ASSUMING points in projection'))}
  sp.xydat <- sf::st_as_sf(xydat,coords=c(xf,yf),crs = .points.crs)
  sf::sf_use_s2(F)
  overlay.sp.xydat <- suppressMessages (as.data.frame (sf::st_join(sp.xydat,.countries.shapefile)))
  country_ext = as.character (overlay.sp.xydat[, ctryNameField])
  country_ext
  }
  

