### Functions slightly modified from other packages for consistency ###
# ne_download_occTest ====
#' @title Download quitely from neearth data (modified from original package)
#' @description download from Natural Earth the country checklist
#' @details
#' @param scale Data.frame of species occurrences
#' @param type the field in the dataframe containing the x coordinates
#' @param category the field in the dataframe containing the y coordinates
#' @param desetdir proj4string argument for dataframe
#' @param load Raster of human influence index
#' @return spatial object or raster
#' @keywords internal
#' @author 
#' @note
#' @seealso rnaturalearth
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
ne_download_occTest = function (scale = 110, type = "countries", 
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



# cc_urb_occTest ====
#' @title Flag values in urban areas (based on coordinate cleaner but modified for occTest to deliver it quitely)
#'
#' @description flags values for urban areas
#' @details
#' @param x Data.frame of species occurrences
#' @param lon character. Column name in x with decimal longitude values
#' @param lat character. Column name in x with decimal latitude values
#' @param a SpatialPolygonsDataFrame. Providing the geographic gazetteer with the urban areas. See details. By default rnaturalearth::ne_download(scale = 'medium', type = 'urban_areas'). Can be any SpatialPolygonsDataframe, but the structure must be identical to rnaturalearth::ne_download(). 
#' @param value character string. Defining the output value. See value.
#' @param verbose logical. If TRUE reports the name of the test and the number of records flagged.
#' @return 
#' @keywords internal
#' @author 
#' @note
#' @seealso CoordinateCleaner::cc_urb
#' @references CoordinateCleaner package
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' 
cc_urb_occTest <-  function (x, lon = "decimallongitude", lat = "decimallatitude", 
                        ref = NULL, value = "clean", verbose = F,outdir=output.dir) 
{
  match.arg(value, choices = c("clean", "flagged"))
  if (verbose) {
    message("Testing urban areas")
  }
  if (is.null(ref)) {
    #message("Downloading urban areas via rnaturalearth")
    ref <- try(suppressWarnings(occTest:::ne_download_occTest(scale = "medium", 
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
    newoutdir = paste0(outdir,'/spatialData')
     dir.create(newoutdir,showWarnings = F,recursive = T)
     rgdal:::writeOGR(ref,dsn = newoutdir,layer = 'NE_urbanareas',driver = "ESRI Shapefile")
     
  }
  else {
    if (!any(is(ref) == "Spatial")) {
      ref <- as(ref, "Spatial")
    }
    #this line is in the original code of coordinate cleaner but I do not know why, I guess it comes from pkg reproj
    #ref <- reproj(ref)
  }
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  dat <- sp::SpatialPoints(x[, c(lon, lat)], proj4string = sp::CRS(wgs84))
  limits <- raster::extent(dat) + 1
  ref <- raster::crop(ref, limits)
  sp::proj4string(ref) <- wgs84
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


# cc_outl_occTest ====
#' @title Identify geographic outliers based on methods from CoordinateCleaner package
#' @description own version of coordinate cleaner geographic outliers
#' @details
#' @keywords internal
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }

cc_outl_occTest <- function (x, lon = "decimallongitude", lat = "decimallatitude", 
                       species = "species", method = "quantile", mltpl = 5, 
                       tdi = 1000, value = "clean", sampling_thresh = 0, verbose = TRUE, 
                       min_occs = 7, thinning = FALSE, thinning_res = 0.5) 
{
  match.arg(value, choices = c("clean", "flagged", 
                               "ids"))
  match.arg(method, choices = c("distance", "quantile", 
                                "mad"))
  if (verbose) {
    message("Testing geographic outliers")
  }
  
  splist <- split(x, f = as.character(x[[species]]))
  test <- lapply(splist, function(k) {
    duplicated(k[, c(species, lon, lat)])
  })
  test <- lapply(test, "!")
  test <- as.vector(unlist(lapply(test, "sum")))
  splist <- splist[test >= min_occs]
  if (any(test < min_occs)) {
    warning(sprintf("Species with less than %o unique records will not be tested.", 
                    min_occs))
  }
  if (any(test >= 10000) | thinning) {
    warning("Using raster approximation.")
    ras <- ras_create(x = x, lat = lat, lon = lon, thinning_res = thinning_res)
  }
  flags <- lapply(splist, function(k) {
    if (nrow(k) <= 10000 & !thinning) {
      dist <- geosphere::distm(k[, c(lon, lat)], fun = geosphere::distHaversine)/1000
      dist[dist == 0] <- NA
    }
    else {
      if (thinning) {
        dist_obj <- ras_dist(x = k, lat = lat, lon = lon, 
                             ras = ras, weights = FALSE)
        pts <- dist_obj$pts
        dist <- dist_obj$dist
      }
      else {
        dist_obj <- ras_dist(x = k, lat = lat, lon = lon, 
                             ras = ras, weights = TRUE)
        pts <- dist_obj$pts
        dist <- dist_obj$dist
        wm <- dist_obj$wm
      }
    }
    if (method == "distance") {
      if (sampling_thresh > 0) {
        stop("Sampling correction impossible for method 'distance'")
      }
      mins <- apply(dist, 1, min, na.rm = TRUE)
      out <- which(mins > tdi)
    }
    else {
      if (nrow(k) >= 10000 & !thinning) {
        mins <- apply(dist, 1, sum, na.rm = TRUE)/rowSums(wm, 
                                                          na.rm = TRUE)
      }
      else {
        mins <- apply(dist, 1, mean, na.rm = TRUE)
      }
    }
    if (method == "quantile") {
      quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE)
      out <- which(mins > quo[2] + stats::IQR(mins) * mltpl)
    }
    if (method == "mad") {
      quo <- stats::median(mins, na.rm = TRUE)
      tester <- stats::mad(mins, na.rm = TRUE)
      out <- which(mins > quo + tester * mltpl)
    }
    if (nrow(k) > 10000 | thinning) {
      if (length(out) == 0) {
        ret <- NA
      }
      if (length(out) > 0) {
        ret <- rownames(k)[which(pts %in% gsub("X", 
                                               "", names(out)))]
      }
    }
    else {
      if (length(out) == 0) {
        ret <- NA
      }
      if (length(out) > 0) {
        ret <- rownames(k)[out]
      }
    }
    return(ret)
  })
  flags <- as.numeric(as.vector(unlist(flags)))
  flags <- flags[!is.na(flags)]
  if (sampling_thresh > 0 & length(flags)>0) {
    pts <- sp::SpatialPoints(x[flags, c(lon, lat)])
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
      stop("package 'rnaturalearth' not found. Needed for sampling_cor = TRUE", 
           call. = FALSE)
    }
    if (!requireNamespace("rgbif", quietly = TRUE)) {
      stop("package 'rgbif' not found. Needed for sampling_cor = TRUE", 
           call. = FALSE)
    }
    if (class(try(rgbif::occ_count(country = "DEU"))) == 
        "try-error") {
      warnings("Could not retrive records number from GBIF, skipping sampling correction")
    }
    else {
      ref <- rnaturalearth::ne_countries(scale = "medium")
      sp::proj4string(ref) <- ""
      #change from AZizka: changed to iso_a2 because in vapply rgbiff:occ count looks for 2 digit country codes
      area <- data.frame(country = ref@data$iso_a2, area = geosphere::areaPolygon(ref))
      area <- area[!is.na(area$area), ]
      area <- area[!is.na(area$country), ]
      nrec <- vapply(area$country, FUN = function(k) {
        rgbif::occ_count(country = k)
      }, FUN.VALUE = 1)
      nrec <- data.frame(country = area$country, recs = unlist(nrec), 
                         row.names = NULL)
      nrec_norm <- dplyr::left_join(nrec, area, by = "country")
      nrec_norm$norm <- log(nrec_norm$recs/(nrec_norm$area/1e+06/100))
      
      #change from AZizka: commented out bc it introduces warnings that may freak out users (e.g. the etent of the coordinates may go beyond the world)
      #ref <- raster::crop(ref, raster::extent(pts) + 1)
      
      country <- sp::over(x = pts, y = ref)[, "iso_a3"]
      thresh <- stats::quantile(nrec_norm$norm, probs = sampling_thresh)
      s_flagged <- nrec_norm$norm[match(country, nrec_norm$country)]
      s_flagged <- s_flagged > thresh
      s_flagged[is.na(s_flagged)] <- FALSE
      flags <- flags[s_flagged]
    }
  }
  out <- rep(TRUE, nrow(x))
  out[flags] <- FALSE
  if (verbose) {
    if (value == "ids") {
      if (value == "clean") {
        message(sprintf("Removed %s records.", 
                        length(flags)))
      }
      else {
        message(sprintf("Flagged %s records.", 
                        length(flags)))
      }
    }
    else {
      if (value == "clean") {
        message(sprintf("Removed %s records.", 
                        sum(!out)))
      }
      else {
        message(sprintf("Flagged %s records.", 
                        sum(!out)))
      }
    }
  }
  switch(value, clean = return(x[out, ]), flagged = return(out), 
         ids = return(flags))
}

