### Functions slightly modified from other packages for consistency with occTest outputs###

# .ne_download_occTest ====
#' @title Download quitely from neearth data (modified from original package)
#' @description download from Natural Earth the country checklist
#' @details
#' @author Andy South (southandy@@gmail.com) rnaturalearth; Josep M Serra Diaz (modifs)
#' @param scale Data.frame of species occurrences
#' @param type the field in the dataframe containing the x coordinates
#' @param category the field in the dataframe containing the y coordinates
#' @param desetdir proj4string argument for dataframe
#' @param load Raster of human influence index
#' @return spatial object or raster
#' @keywords internal
#' @seealso \[linkrnaturalearth]{ne_download}
#' @import rnaturalearth

.ne_download_occTest = function (scale = 110, type = "countries", 
                                    category = c("cultural", "physical", "raster"), 
                                    destdir = tempdir(), 
                                    load = TRUE, 
                                    returnclass = c("sp", "sf")) {
  category <- match.arg(category)
  returnclass <- match.arg(returnclass)
  file_name <- rnaturalearth::ne_file_name(scale = scale, type = type, category = category, 
                                            full_url = FALSE)
  address <- rnaturalearth::ne_file_name(scale = scale, type = type, category = category, 
                                          full_url = TRUE)
  utils::download.file(file.path(address), zip_file <- tempfile(),quiet = TRUE)
  utils::unzip(zip_file, exdir = destdir)
  if (load & category == "raster") {
    rst <- raster::raster(file.path(destdir, file_name, paste0(file_name, 
                                                                ".tif")))
    return(rst)
  }
  else if (load) {
    sp_object <- rgdal::readOGR(destdir, file_name, encoding = "UTF-8", 
                                 stringsAsFactors = FALSE, use_iconv = TRUE,verbose = FALSE)
    sp_object@data[sp_object@data == "-99" | sp_object@data == 
                     "-099"] <- NA
    return(ne_as_sf(sp_object, returnclass))
  }
  else {
    return(file_name)
  }
}



# .cc_urb_occTest ====
#' @title Flag values in urban areas (based on coordinate cleaner but modified for occTest to deliver it quitely)
#' @description flags values for urban areas
#' @param x Data.frame of species occurrences
#' @param lon character. Column name in x with decimal longitude values
#' @param lat character. Column name in x with decimal latitude values
#' @param a SpatialPolygonsDataFrame. Providing the geographic gazetteer with the urban areas. See details. By default rnaturalearth::ne_download(scale = 'medium', type = 'urban_areas'). Can be any SpatialPolygonsDataframe, but the structure must be identical to rnaturalearth::ne_download(). 
#' @param value character string. Defining the output value. See value.
#' @param verbose logical. If TRUE reports the name of the test and the number of records flagged
#' @param outdir output directory
#' @return a clean data.frame 
#' @keywords internal
#' @author A Zizka (original function) Josep M Serra-Diaz (adaptation to occTest pep.serradiaz@@agroparistech.fr)
#' @seealso \link[CoordinateCleaner]{cc_urb}
#' @references CoordinateCleaner package
#' @importFrom methods as is slot
.cc_urb_occTest <-  function (x, lon = "decimallongitude", lat = "decimallatitude", 
                        ref = NULL, value = "clean", verbose = FALSE,outdir) 
{
  match.arg(value, choices = c("clean", "flagged"))
  if (verbose) {
    message("Testing urban areas")
  }
  if (is.null(ref)) {
    #message("Downloading urban areas via rnaturalearth")
    ref <- try(suppressWarnings(.ne_download_occTest(scale = "medium", 
                                                 type = "urban_areas")), silent = TRUE)
    if (inherits(ref,"try-error") ) {
      warning(sprintf("Gazetteer for urban areas not found at\n%s", 
                      rnaturalearth::ne_file_name(scale = "medium", 
                                                  type = "urban_areas", full_url = TRUE)))
      warning("Skipping urban test")
      switch(value, clean = return(x), flagged = return(rep(NA, 
                                                            nrow(x))))
    }
    sp::proj4string(ref) <- ""
    newoutdir = paste0(outdir,'/spatialData')
     dir.create(newoutdir,showWarnings = FALSE,recursive = TRUE)
     rgdal::writeOGR(ref,dsn = newoutdir,layer = 'NE_urbanareas',driver = "ESRI Shapefile",overwrite_layer=TRUE)
     
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
  sp::proj4string(ref) <- wgs84
  ref <- raster::crop(ref, limits)
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


# .cc_outl_occTest ====
#' @title Identify geographic outliers based on methods from CoordinateCleaner package
#' @description own version of coordinate cleaner geographic outliers
#' @param x Data.frame of species occurrences
#' @param lon character. Column name in x with decimal longitude values
#' @param lat character. Column name in x with decimal latitude values
#' @param species  character. Species name. 
#' @param method character. Quantile by default
#' @param mltpl integer. Multiplication factor. Default to 5
#' @param tdi numeric. Default to 1000
#' @param value character. Type of output. Default to 'clean'
#' @param sampling_thresh numeric. Sampling threshold. Defaults to 0
#' @param verbose logical. Should output messages be printed. Default to TRUE
#' @param min_occs integer. Minimum number of occurrences. Defaults to 7
#' @param thinning logical. Should thinning be performed? Defaults to FALSE
#' @param thinning_res double. Thinnning resolution. Defaults to 0.5
#' @seealse \link[CoordinateCleaner]{cc_outl}
#' @keywords internal
#' @author A Zizka (original function) Josep M Serra-Diaz (adaptation to occTest pep.serradiaz@@agroparistech.fr)
#' @return a clean data.frame 
#' @import CoordinateCleaner
.cc_outl_occTest <- function (x, lon = "decimallongitude", lat = "decimallatitude", 
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
      out <- which(mins > quo[2] + stats::IQR(mins,na.rm = TRUE) * mltpl)
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
    if (inherits(try(rgbif::occ_count(country = "DEU")),"try-error")) {
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


# cc_round_occTest ====
#' @title Flag records with regular pattern interval
#' @description own version of coordinate cleaner cc_round
#' @author A Zizka (original author) Josep M Serra-Diaz (adapted from CoordinateCleaner)
#' @param x Data.frame of species occurrences
#' @param lon character. Column name in x with decimal longitude values
#' @param lat character. Column name in x with decimal latitude values
#' @param ds character. Column name in x with dataset name of the record
#' @param T1 numeric. Defaults to 7
#' @param reg_out_thresh numeric. Defaults to 7
#' @param reg_dist_min numeric. Defaults to 7
#' @param reg_dist_max numeric. Defaults to 7
#' @param min_unique_ds_size numeric. Defaults to 7
#' @param graphs logical. Defaults to FALSE.
#' @param test character. Defaults to 'both'
#' @param value character. Defaults to flagged
#' @param verbose logical. Defaults to TRUE.
#' @seealso \link[CoordinateCleaner]{CoordinateCleaner-package}
#' @notes Turned off by default as it disappeared from CoordinateCleaner pkg
#' @return a clean data.frame 
#' @import CoordinateCleaner
#' @importFrom graphics title
.cc_round_occTest <-  function (x, lon = "decimallongitude", lat = "decimallatitude", 
                           ds = "dataset", T1 = 7, reg_out_thresh = 2, reg_dist_min = 0.1, 
                           reg_dist_max = 2, min_unique_ds_size = 4, graphs = FALSE, 
                           test = "both", value = "flagged", verbose = TRUE) 
{
  window_size <- 10
  detection_rounding <- 2
  detection_threshold <- 6
  digit_round <- 0
  nc <- 3000
  rarefy <- FALSE
  match.arg(value, choices = c("flagged", "clean", "dataset","flagged2"))
  if (verbose) {
    message("Testing for rasterized collection")
  }
  if (length(unique(x[[ds]])) > 1) {
    dat <- split(x, f = x[[ds]])
    out <- lapply(dat, function(k) {


      tester <- k[stats::complete.cases(k[, c(lon, lat)]), ]
      if (nrow(tester[!duplicated(tester[, c(lon, lat)]), 
      ]) < min_unique_ds_size) {
        warning(paste0 (unique(k[[ds]])," :Dataset smaller than minimum test size"))
        # out <- data.frame(dataset = unique(x[[ds]]), 
        #                   n.outliers = NA, n.regular.outliers = NA, 
        #                   regular.distance = NA, summary = NA)
        out <- data.frame(dataset = unique(k[[ds]]),
                          n.outliers = NA, n.regular.outliers = NA,
                          regular.distance = NA, summary = NA)

      }
      else {
        if (test == "lon") {
          gvec <- try(.CalcACT(data = k[[lon]], digit_round = digit_round, 
                           nc = nc, graphs = graphs, graph_title = unique(k[[ds]])),silent=TRUE)
          if (class(gvec) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl <- try(.OutDetect(gvec, T1 = T1, window_size = window_size, 
                               detection_rounding = detection_rounding, 
                               detection_threshold = detection_threshold, 
                               graphs = graphs),silent=TRUE)
          
          if (class(n_outl) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl$flag <- !all(n_outl$n.outliers > 0, 
                              n_outl$regular.distance >= reg_dist_min, 
                              n_outl$regular.distance <= reg_dist_max, 
                              n_outl$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl$flag, 
                        sep = " - "))
          }
          n_outl <- data.frame(unique(k[[ds]]), n_outl)
          names(n_outl) <- c("dataset", "lon.n.outliers", 
                             "lon.n.regular.distance", "lon.regular.distance", 
                             "summary")
        }
        if (test == "lat") {
          gvec <- try({.CalcACT(data = k[[lat]], digit_round = digit_round, 
                           nc = nc, graphs = graphs, graph_title = unique(k[[ds]]))}, 
                      silent=TRUE)
          
          if (class(gvec) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl <- try({.OutDetect(gvec, T1 = T1, window_size = window_size, 
                               detection_rounding = detection_rounding, 
                               detection_threshold = detection_threshold, 
                               graphs = graphs)},silent = TRUE)
          
          if (class(n_outl) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl$flag <- !all(n_outl$n.outliers > 0, 
                              n_outl$regular.distance >= reg_dist_min, 
                              n_outl$regular.distance <= reg_dist_max, 
                              n_outl$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl$flag, 
                        sep = " - "))
          }
          n_outl <- data.frame(unique(k[[ds]]), n_outl)
          names(n_outl) <- c("dataset", "lat.n.outliers", 
                             "lat.n.regular.distance", "lat.regular.distance", 
                             "summary")
        }
        
        if (test== 'both') {
          gvec1 <- try(.CalcACT(data = k[[lon]], digit_round = digit_round, 
                                                   nc = nc, graphs = graphs, graph_title = unique(k[[ds]])),silent=TRUE)
          if (class(gvec1) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl_lon <- try(.OutDetect(gvec1, T1 = T1, window_size = window_size, 
                                                       detection_rounding = detection_rounding, 
                                                       detection_threshold = detection_threshold, 
                                                       graphs = graphs),silent=TRUE)
          
          if (class(n_outl_lon) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl_lon$flag <- !all(n_outl_lon$n.outliers > 0, 
                                  n_outl_lon$regular.distance >= reg_dist_min, 
                                  n_outl_lon$regular.distance <= reg_dist_max, 
                                  n_outl_lon$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl_lon$flag, 
                        sep = " - "))
          }
          
          #latitude
          gvec2 <- try({.CalcACT(data = k[[lat]], digit_round = digit_round, 
                                                    nc = nc, graphs = graphs, graph_title = unique(k[[ds]]))}, 
                      silent=TRUE)
          
          if (class(gvec2) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl_lat <- try({.OutDetect(gvec2, T1 = T1, window_size = window_size, 
                                                        detection_rounding = detection_rounding, 
                                                        detection_threshold = detection_threshold, 
                                                        graphs = graphs)},silent = TRUE)
          
          if (class(n_outl_lat) %in% c('error','try-error')) {
            out <- data.frame(dataset = unique(k[[ds]]),
                              n.outliers = NA, n.regular.outliers = NA,
                              regular.distance = NA, summary = NA)
            return (out)
          }
          
          n_outl_lat$flag <- !all(n_outl_lat$n.outliers > 0, 
                                  n_outl_lat$regular.distance >= reg_dist_min, 
                                  n_outl_lat$regular.distance <= reg_dist_max, 
                                  n_outl_lat$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl_lat$flag, 
                        sep = " - "))
          }
          
          n_outl <- cbind(unique(k[[ds]]), n_outl_lon, 
                          n_outl_lat)
          names(n_outl) <- c("dataset", "lon.n.outliers", 
                             "lon.n.regular.outliers", "lon.regular.distance", 
                             "lon.flag", "lat.n.outliers", "lat.n.regular.outliers", 
                             "lat.regular.distance", "lat.flag")
          n_outl$summary <- n_outl$lon.flag | n_outl$lat.flag
        }
        if (test == "bothOld") {
          gvec1 <- .CalcACT(data = k[[lon]], digit_round = digit_round, 
                            nc = nc, graphs = graphs, graph_title = unique(k[[ds]]))
          n_outl_lon <- .OutDetect(gvec1, T1 = T1, window_size = window_size, 
                                   detection_rounding = detection_rounding, 
                                   detection_threshold = detection_threshold, 
                                   graphs = graphs)
          n_outl_lon$flag <- !all(n_outl_lon$n.outliers > 
                                    0, n_outl_lon$regular.distance >= reg_dist_min, 
                                  n_outl_lon$regular.distance <= reg_dist_max, 
                                  n_outl_lon$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl_lon$flag, 
                        sep = " - "))
          }
          gvec2 <- .CalcACT(data = k[[lat]], digit_round = digit_round, 
                            nc = nc, graphs = graphs, graph_title = unique(k[[ds]]))
          n_outl_lat <- .OutDetect(gvec2, T1 = T1, window_size = window_size, 
                                   detection_rounding = detection_rounding, 
                                   detection_threshold = detection_threshold, 
                                   graphs = graphs)
          n_outl_lat$flag <- !all(n_outl_lat$n.outliers > 
                                    0, n_outl_lat$regular.distance >= reg_dist_min, 
                                  n_outl_lat$regular.distance <= reg_dist_max, 
                                  n_outl_lat$n.regular.outliers >= reg_out_thresh)
          if (graphs) {
            title(paste(unique(k[[ds]]), n_outl_lat$flag, 
                        sep = " - "))
          }
          n_outl <- cbind(unique(k[[ds]]), n_outl_lon, 
                          n_outl_lat)
          names(n_outl) <- c("dataset", "lon.n.outliers", 
                             "lon.n.regular.outliers", "lon.regular.distance", 
                             "lon.flag", "lat.n.outliers", "lat.n.regular.outliers", 
                             "lat.regular.distance", "lat.flag")
          n_outl$summary <- n_outl$lon.flag | n_outl$lat.flag
        }
        return(n_outl)
      }
    })
    #out <- do.call("rbind.data.frame", out)

    out <- do.call("bind_rows", out)
  }
  else {
    if (nrow(x[!duplicated(x[, c(lon, lat)]), ]) < min_unique_ds_size) {
      warning("Dataset smaller than minimum test size")
      out <- data.frame(dataset = unique(x[[ds]]), n.outliers = NA, 
                        n.regular.outliers = NA, regular.distance = NA, 
                        summary = NA)
    }
    else {
      if (test == "lon") {
        gvec <- try({.CalcACT(data = x[[lon]], digit_round = digit_round, 
                         nc = nc, graphs = graphs, graph_title = unique(x[[ds]]))},silent=TRUE)
        if (class(gvec) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
        }
        
        n_outl <- .OutDetect(gvec, T1 = T1, window_size = window_size, 
                             detection_rounding = detection_rounding, detection_threshold = detection_threshold, 
                             graphs = graphs)
        
        if (class(n_outl) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )        
          }
        
        n_outl$flag <- !all(n_outl$n.outliers > 0, n_outl$regular.distance >= 
                              reg_dist_min, n_outl$regular.distance <= reg_dist_max, 
                            n_outl$n.regular.outliers >= reg_out_thresh)
        if (graphs) {
          title(paste(unique(x[[ds]]), n_outl$flag, 
                      sep = " - "))
        }
        n_outl <- data.frame(unique(x[[ds]]), n_outl)
        names(n_outl) <- c("dataset", "lon.n.outliers", 
                           "lon.n.regular.distance", "lon.regular.distance", 
                           "summary")
      }
      if (test == "lat") {
        gvec <- try (.CalcACT(data = x[[lat]], digit_round = digit_round, 
                         nc = nc, graphs = graphs, graph_title = unique(x[[ds]])),silent=TRUE)
        if (class(gvec) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
          }
        n_outl <- try (.OutDetect(gvec, T1 = T1, window_size = window_size, 
                             detection_rounding = detection_rounding, detection_threshold = detection_threshold, 
                             graphs = graphs), silent=TRUE)
        if (class(n_outl) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )        
          }
        
        n_outl$flag <- !all(n_outl$n.outliers > 0, n_outl$regular.distance >= 
                              reg_dist_min, n_outl$regular.distance <= reg_dist_max, 
                            n_outl$n.regular.outliers >= reg_out_thresh)
        if (graphs) {
          title(paste(unique(x[[ds]]), n_outl$flag, 
                      sep = " - "))
        }
        n_outl <- data.frame(unique(x[[ds]]), n_outl)
        names(n_outl) <- c("dataset", "lat.n.outliers", 
                           "lat.n.regular.distance", "lat.regular.distance", 
                           "summary")
      }
      if (test == "both") {
        gvec1 <- try(.CalcACT(data = x[[lon]], digit_round = digit_round, 
                          nc = nc, graphs = graphs, graph_title = unique(x[[ds]])), silent=TRUE)
        
        if (class(gvec1) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
          }
        
        n_outl_lon <- try (.OutDetect(gvec1, T1 = T1, window_size = window_size, 
                                 detection_rounding = detection_rounding, detection_threshold = detection_threshold, 
                                 graphs = graphs), silent=TRUE)
        if (class(n_outl_lon) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
          }
        
        
        n_outl_lon$flag <- !all(n_outl_lon$n.outliers > 
                                  0, n_outl_lon$regular.distance >= reg_dist_min, 
                                n_outl_lon$regular.distance <= reg_dist_max, 
                                n_outl_lon$n.regular.outliers >= reg_out_thresh)
        if (graphs) {
          title(paste(unique(x[[ds]]), n_outl_lon$flag, 
                      sep = " - "))
        }
        gvec2 <- try (.CalcACT(data = x[[lat]], digit_round = digit_round, 
                          nc = nc, graphs = graphs, graph_title = unique(x[[ds]])), silent=TRUE)
        if (class(gvec2) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
          }
        n_outl_lat <- try (.OutDetect(gvec2, T1 = T1, window_size = window_size, 
                                 detection_rounding = detection_rounding, detection_threshold = detection_threshold, 
                                 graphs = graphs), silent = TRUE)
        
        if (class(n_outl_lat) %in% c('error','try-error')) {
          out <- data.frame(dataset = rep(unique(x[[ds]]),length.out=nrow(x)), 
                            n.outliers = rep(NA,length.out=nrow(x)), 
                            n.regular.outliers = rep(NA,length.out=nrow(x)), 
                            regular.distance = rep(NA,length.out=nrow(x)), 
                            summary = rep(NA,length.out=nrow(x)))
          
          switch(value, dataset = return(out), clean = return({
            test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
            ]
            if (length(test) > 0) {
              test
            } else {
              NULL
            }
          }), 
          flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
          flagged2 = return(out$summary)
          )
          }
        
        n_outl_lat$flag <- !all(n_outl_lat$n.outliers > 
                                  0, n_outl_lat$regular.distance >= reg_dist_min, 
                                n_outl_lat$regular.distance <= reg_dist_max, 
                                n_outl_lat$n.regular.outliers >= reg_out_thresh)
        if (graphs) {
          title(paste(unique(x[[ds]]), n_outl_lat$flag, 
                      sep = " - "))
        }
        n_outl <- data.frame(unique(x[[ds]]), n_outl_lon, 
                             n_outl_lat)
        names(n_outl) <- c("dataset", "lon.n.outliers", 
                           "lon.n.regular.distance", "lon.regular.distance", 
                           "lon.flag", "lat.n.outliers", "lat.n.regular.distance", 
                           "lat.regular.distance", "lat.flag")
        n_outl$summary <- n_outl$lon.flag | n_outl$lat.flag
      }
      out <- n_outl
    }
  }
  out2 = merge (x,out[c('dataset','summary')], by.x= ds,by.y='dataset',all.x=TRUE)
  switch(value, dataset = return(out), clean = return({
    test <- x[x[[ds]] %in% out[out$summary, "dataset"], 
    ]
    if (length(test) > 0) {
      test
    } else {
      NULL
    }
  }), 
  flagged = return(x[[ds]] %in% out[out$summary, "dataset"]),
  flagged2 = return(out2$summary)
  )
}


## .cd_ddmm occTest ====
#' @title Flag records with regular pattern interval
#' @description own version of coordinate cleaner cc_round
#' @details
#' @keywords internal
#' @author A Zizka (original author) Josep M Serra-Diaz (adapted from CoordinateCleaner)
#' @return a clean data.frame 

.cd_ddmm_occTest <- function (x, lon = "decimallongitude", lat = "decimallatitude", 
          ds = "dataset", pvalue = 0.025, diff = 1, mat_size = 1000, 
          min_span = 2, value = "clean", verbose = TRUE, diagnostic = FALSE) 
{
  match.arg(value, choices = c("clean", "flagged", "dataset"))
  if (verbose) {
    message("Testing for dd.mm to dd.dd conversion errors")
  }
  if (sum(!stats::complete.cases(x[, c(lon, lat, ds)])) > 0) {
    warning(sprintf("ignored %s cases with incomplete data", 
                    sum(!stats::complete.cases(x))))
  }
  dat <- x[stats::complete.cases(x[, c(lon, lat, ds)]), ]
  if (nrow(dat) == 0) {
    stop("no complete cases found")
  }
  dat$lon.test <- abs(dat[[lon]]) - floor(abs(dat[[lon]]))
  dat$lat.test <- abs(dat[[lat]]) - floor(abs(dat[[lat]]))
  test <- split(dat, f = dat[[ds]])
  out <- lapply(test, function(k) {
    dat_unique <- k[!duplicated(k[, c(lon, lat, ds)]), ]
    lon_span <- abs(max(dat_unique[, lon], na.rm = TRUE) - 
                      min(dat_unique[, lon], na.rm = TRUE))
    lat_span <- abs(max(dat_unique[, lat], na.rm = TRUE) - 
                      min(dat_unique[, lat], na.rm = TRUE))
    if (lon_span >= min_span & lat_span >= min_span) {
      cl <- ceiling(dat_unique[, c("lon.test", "lat.test")] * 
                      mat_size)
      cl$lat.test <- mat_size - cl$lat.test
      mat <- matrix(ncol = mat_size, nrow = mat_size)
      mat[cbind(cl$lat.test, cl$lon.test)] <- 1
      mat[is.na(mat)] <- 0
      dat_t1 <- mat
      P_smaller_than_06 <- floor(0.599 * mat_size) * floor(0.599 * 
                                                             mat_size)/mat_size^2
      x_ind <- (mat_size - floor(0.599 * mat_size)):mat_size
      y_ind <- 1:floor(0.599 * mat_size)
      subt <- dat_t1[x_ind, y_ind]
      p06 <- sum(subt >= 1)
      pAll <- sum(dat_t1 >= 1)
      #print(unique (k[,ds]))
      if (pAll <1)  {outp <- rep(NA, 3) ; return(outp)}
      B <- stats::binom.test(p06, pAll, p = P_smaller_than_06, 
                             alternative = c("greater"))
      v1 <- B$p.value
      v2 <- (B$estimate - P_smaller_than_06)/P_smaller_than_06
      if (v1 < pvalue & v2 > diff) {
        flag_t1 <- FALSE
      }
      else {
        flag_t1 <- TRUE
      }
      outp <- c(round(v1, 4), round(v2, 3), flag_t1)
      if (diagnostic) {
        plo <- raster::raster(dat_t1)
        raster::plot(plo)
        title(as.logical(flag_t1))
      }
    }
    else {
      outp <- rep(NA, 3)
      warning("Geographic span too small, check 'min_span'")
    }
    return(outp)
  })
  out_ds <- do.call("rbind.data.frame", out)
  rownames(out_ds) <- names(out)
  names(out_ds) <- c("binomial.pvalue", "perc.difference", 
                     "pass")
  flags <- x[[ds]] %in% c(rownames(out_ds[out_ds$pass == 1, 
  ]), rownames(out_ds[is.na(out_ds$pass), ]))
  out_ds$pass <- as.logical(out_ds$pass)
  if (verbose) {
    message(sprintf("Flagged %s records", sum(!flags)))
  }
  switch(value, dataset = return(out_ds), clean = return(x[flags, 
  ]), flagged = return(flags))
}





