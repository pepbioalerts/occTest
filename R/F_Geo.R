#' Identifies country ISO3 based on locations
#'
#' extracts data from database 
#' @param points A dataframe with x and y coordinates
#' @return Factor with ISO3 codes for countries 
#' @family Geo
#' @examples \dontrun{
#'
#' }
#' 
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
.coords2country = function(points)
{  
  countriesSP <- rworldmap:::getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = sp:::SpatialPoints(points, proj4string=sp:::CRS(sp:::proj4string(countriesSP)))  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = sp:::over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  #indices$ADMIN  
  as.character(indices$ISO3) # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}



#' Extracts a dataframe based on dataframe of xy locations
#'
#' Based on getData from raster
#' @param xydat A dataframe with x and y coordinates
#' @param download Default to T. Whether the data should be downloaded 
#' @param path where the downloads should go. Default to the currendt directory
#' @param verbose if you want to print messages of progress or warnings
#' @param 
#' @return Factor with ISO3 codes for countries 
#' @family Geo
#' @examples \dontrun{
#'
#' }
#' 

.getSRTM = function (xydat,download=T,path=getwd(), verbose=F) {
  #get file names of the tiles in srtm
  file.tiles = sapply ( 1:nrow (xydat), function (i) {
    lon <- as.numeric (xydat[i,1])
    lat <- as.numeric (xydat[i,1])
    stopifnot(lon >= -180 & lon <= 180)
    stopifnot(lat >= -60 & lat <= 60)
    rs <- raster:::raster(nrows = 24, ncols = 72, xmn = -180, xmx = 180, 
                          ymn = -60, ymx = 60)
    rowTile <- raster:::rowFromY(rs, lat)
    colTile <- raster:::colFromX(rs, lon)
    if (rowTile < 10) {
      rowTile <- paste("0", rowTile, sep = "")
    }
    if (colTile < 10) {
      colTile <- paste("0", colTile, sep = "")
    }
    f <- paste0("srtm_", colTile, "_", rowTile, ".zip")
    return (f)
    
  })
  file.tiles <- unique (file.tiles)
  
  #download ouptus
  ras.list = lapply (file.tiles, function (f){
    baseurl <- "http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/"
    zipfilename <- file.path(path, f)
    tiffilename <- file.path(path, gsub(".zip$", ".tif", f))
    if (!file.exists(tiffilename)) {
      if (!file.exists(zipfilename)) {
        if (download) {
          theurl <- paste0(baseurl, f)
          test <- try(raster:::.download(theurl, zipfilename), silent = TRUE)
          if (class(test) == "try-error") {
            stop("cannot download the file")
          }
        }
        else {
          message("file not available locally, use download=TRUE")
        }
      }
      if (file.exists(zipfilename)) {
        utils::unzip(zipfilename, exdir = dirname(zipfilename))
        file.remove(zipfilename)
      }
    }
    if (file.exists(tiffilename)) {
      rs <- raster(tiffilename)
      raster:::projection(rs) <- "+proj=longlat +datum=WGS84"
      return(rs)
    }
    else {
      stop("file not found")
    }
  })
  
  #merge rasters
  if (length (ras.list) >1) {srtm.ras <- do.call(raster::merge, ras.list)}
  if (length (ras.list) == 1) {srtm.ras <- ras.list[[1]]}
  
  #gimme the raster back
  srtm.ras
  
}
