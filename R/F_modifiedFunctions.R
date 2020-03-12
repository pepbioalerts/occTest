#' @title Download (QUIETELY!) from neearth data (modified from original package)
#'
#' @description download from Natural Earth
#' @details
#' @param scale Data.frame of species occurrences
#' @param type the field in the dataframe containing the x cordinates
#' @param category the field in the dataframe containing the y cordinates
#' @param desetdir proj4string argument for dataframe
#' @param load Raster of human influence index
#' @return spatial object or raster
#' @keywords internal
#' @author author of the package neearth
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }


ne_download_occProfileR = function (scale = 110, type = "countries", 
                                    category = c("cultural", "physical", "raster"), 
                                    destdir = tempdir(), 
                                    load = TRUE, 
                                    returnclass = c("sp", "sf")) {
  category <- match.arg(category)
  returnclass <- match.arg(returnclass)
  file_name <- rnaturalearth:::ne_file_name(scale = scale, type = type, category = category, 
                                            full_url = FALSE)
  address <- rnaturalearth:::ne_file_name(scale = scale, type = type, category = category, 
                                          full_url = TRUE)
  utils:::download.file(file.path(address), zip_file <- tempfile(),quiet = T)
  utils:::unzip(zip_file, exdir = destdir)
  if (load & category == "raster") {
    rst <- raster:::raster(file.path(destdir, file_name, paste0(file_name, 
                                                                ".tif")))
    return(rst)
  }
  else if (load) {
    sp_object <- rgdal:::readOGR(destdir, file_name, encoding = "UTF-8", 
                                 stringsAsFactors = FALSE, use_iconv = TRUE,verbose = F)
    sp_object@data[sp_object@data == "-99" | sp_object@data == 
                     "-099"] <- NA
    return(rnaturalearth:::ne_as_sf(sp_object, returnclass))
  }
  else {
    return(file_name)
  }
}




#' @title Flag values in urban areas (based on coordinate cleaner but modified for occProfileR to deliver it quitely)
#'
#' @description flags values for urban areas
#' @details
#' @param x Data.frame of species occurrences
#' @param lon 
#' @param lat 
#' @param ref 
#' @param value 
#' @param verbose 
#' @return 
#' @keywords internal
#' @author author of the coordinatecleaner package (Zizka)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' 
cc_urb_occProfileR <-  function (x, lon = "decimallongitude", lat = "decimallatitude", 
                        ref = NULL, value = "clean", verbose = F) 
{
  match.arg(value, choices = c("clean", "flagged"))
  if (verbose) {
    message("Testing urban areas")
  }
  if (is.null(ref)) {
    #message("Downloading urban areas via rnaturalearth")
    ref <- try(suppressWarnings(occProfileR:::ne_download_occProfileR(scale = "medium", 
                                                 type = "urban_areas")), silent = TRUE)
    if (class(ref) == "try-error") {
      warning(sprintf("Gazetteer for urban areas not found at\n%s", 
                      rnaturalearth::ne_file_name(scale = "medium", 
                                                  type = "urban_areas", full_url = TRUE)))
      warning("Skipping urban test")
      switch(value, clean = return(x), flagged = return(rep(NA, 
                                                            nrow(x))))
    }
    sp::proj4string(ref) <- ""
  }
  else {
    if (!any(is(ref) == "Spatial")) {
      ref <- as(ref, "Spatial")
    }
    ref <- reproj(ref)
  }
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  dat <- sp::SpatialPoints(x[, c(lon, lat)], proj4string = CRS(wgs84))
  limits <- raster::extent(dat) + 1
  ref <- raster::crop(ref, limits)
  proj4string(ref) <- wgs84
  if (is.null(ref)) {
    out <- rep(TRUE, nrow(x))
  }
  else {
    out <- is.na(sp::over(x = dat, y = ref)[, 1])
  }
  if (verbose) {
    if (value == "clean") {
      message(sprintf("Removed %s records.", sum(!out)))
    }
    else {
      message(sprintf("Flagged %s records.", sum(!out)))
    }
  }
  switch(value, clean = return(x[out, ]), flagged = return(out))
}
