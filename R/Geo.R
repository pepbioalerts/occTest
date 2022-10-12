### COLLECTION OF GEOGRAPHIC INFORMATION FUNCTIONS ###

#.coords2country =====
#' @title Extracts country ISO3 based on locations
#' @param xydat A dataframe with x and y coordinates
#' @param .countries.shapefile SpatialPolygonsDataFrame of world countries and their associated ISO3 codes
#' @param .points.proj4string Proj4string for the occurrence data
#' @param ctryNameField character. Column name in .countries shapefile where ISO3 are indicated
#' @param verbose logical. Print messages?
#' @return Factor with ISO3 codes for countries 
#' @family Geo
#' @author Josep M Serra Diaz


.coords2country = function(xydat,
                           .countries.shapefile=NULL,
                           .points.proj4string=NULL,
                           ctryNameField=NULL,verbose=FALSE){  
  if (is.null(.countries.shapefile)) {   .countries.shapefile = rworldmap::getMap(resolution='high') ;  ctryNameField ='ISO3'}
  if (!inherits(.countries.shapefile,'SpatialPolygonsDataFrame')) {stop (".countries shapefile not loaded ")}
  if (is.null(.points.proj4string)) {.points.proj4string <- .countries.shapefile@proj4string; if(verbose) print (paste ('ASSUMING points in projection',.countries.shapefile@proj4string))}
  sp.xydat = sp::SpatialPoints(xydat,proj4string = .points.proj4string)
  overlay.sp.xydat = sp::over(sp.xydat, .countries.shapefile)
  country_ext = as.character (overlay.sp.xydat[, ctryNameField])
  country_ext
  }
  

#.getSRTM =====
#' @title Download SRTM elevation raster
#' @param xydat A dataframe with x and y coordinates
#' @param download Default to TRUE. Whether the data should be downloaded 
#' @param path where the downloads should go. Default to the current directory
#' @param verbose if you want to print messages of progress or warnings
#' @details Basedd on getData from raster
#' @return raster
#' @family Geo
#' @note borrowed from raster package but adapted to work directly within the occTest workflow
#' @seealso  \link[raster]{getData}


.getSRTM = function (xydat,download=TRUE,path=tempdir(), verbose=FALSE) {
  #get file names of the tiles in srtm
  file.tiles = sapply ( 1:nrow (xydat), function (i) {
    lon <- as.numeric (xydat[i,1])
    lat <- as.numeric (xydat[i,1])
    stopifnot(lon >= -180 & lon <= 180)
    stopifnot(lat >= -60 & lat <= 60)
    rs <- raster::raster(nrows = 24, ncols = 72, xmn = -180, xmx = 180, 
                          ymn = -60, ymx = 60)
    rowTile <- raster::rowFromY(rs, lat)
    colTile <- raster::colFromX(rs, lon)
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
  
  #download outputs
  ras.list = lapply (file.tiles, function (f){
    baseurl <- "http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/"
    zipfilename <- file.path(path, f)
    tiffilename <- file.path(path, gsub(".zip$", ".tif", f))
    if (!file.exists(tiffilename)) {
      if (!file.exists(zipfilename)) {
        if (download) {
          theurl <- paste0(baseurl, f)
          test <- try(.download(theurl, zipfilename), silent = TRUE)
          if (inherits(test,"try-error")) {
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
      rs <- raster::raster(tiffilename)
      raster::projection(rs) <- "+proj=longlat +datum=WGS84"
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
