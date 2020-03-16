#' @title POTENTIAL BOTANIC GARDEN LOCALITY FROM NAME
#'
#' @description own version of coordinate cleaner
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
#' @export

Mycc_outl <- function (x, lon = "decimallongitude", lat = "decimallatitude", 
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
