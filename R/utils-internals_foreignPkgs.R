### IMPORTED COPIES OF HIDDEN FUNCTIONS IN OTHER PACKGES (for stability)

# .CalcACT Non-exported function from CoordinateCleaner ====
#' @importFrom graphics hist plot abline segments
#' @importFrom stats cov quantile aggregate
#' @noRd
.CalcACT <- function (data, digit_round = 0, nc = 3000, graphs = TRUE, graph_title = "Title", 
                      rarefy = FALSE) {
  if (rarefy) {
    data <- unique(data)
  }
  data_units <- sort(abs(data))
  if (digit_round > 0) {
    data_units <- data_units - round(floor(data_units/digit_round) * 
                                       digit_round)
  }
  if (graphs) {
    h <- graphics::hist(data_units, nclass = nc, main = graph_title)
  }
  else {
    h <- graphics::hist(data_units, nclass = nc, plot = FALSE)
  }
  f <- h$counts
  max_range <- round(length(f) * 0.9)
  gamma_0 <- stats::cov(f[1:max_range], f[1:max_range])
  gamma_vec <- c()
  for (k in 1:max_range) {
    f_0 <- f[-(1:k)]
    f_k <- f[-((length(f) - k + 1):(length(f)))]
    gamma_vec <- c(gamma_vec, stats::cov(f_0, f_k)/gamma_0)
  }
  coords <- h$mids[1:max_range]
  if (graphs) {
    plot(coords, gamma_vec)
  }
  out <- data.frame(gamma = gamma_vec, coords = coords)
  return(out)
}

# .OutDetect Non-exported function from CoordinateCleaner ====
#' @noRd
.OutDetect <- function (x, T1 = 7, window_size = 10, detection_rounding = 2, 
          detection_threshold = 6, graphs = TRUE) {
  max_range <- nrow(x) - window_size
  out <- matrix(ncol = 2)
  for (k in 1:max_range) {
    sub <- x[k:(k + window_size), ]
    quo <- stats::quantile(sub$gamma, c(0.25, 0.75), na.rm = TRUE)
    outl <- matrix(nrow = nrow(sub), ncol = 2)
    outl[, 1] <- sub$gamma > (quo[2] + stats::IQR(sub$gamma) * 
                                T1)
    outl[, 2] <- sub$coord
    out <- rbind(out, outl)
  }
  out <- stats::aggregate(out[, 1] ~ out[, 2], FUN = "sum")
  names(out) <- c("V1", "V2")
  out <- out[, c(2, 1)]
  out[, 1] <- as.numeric(out[, 1] >= detection_threshold)
  outl <- out[out[, 1] == 1, ]
  if (nrow(outl) == 0) {
    out <- data.frame(n.outliers = 0, n.regular.outliers = 0, 
                      regular.distance = NA)
  }
  else if (nrow(outl) == 1) {
    if (graphs) {
      graphics::abline(v = outl[, 2], col = "green")
    }
    out <- data.frame(n.outliers = 1, n.regular.outliers = 0, 
                      regular.distance = NA)
  }
  else {
    if (graphs) {
      graphics::abline(v = outl[, 2], col = "green")
    }
    dist_m <- round(stats::dist(round(outl[, 2, drop = FALSE], 
                                      detection_rounding), diag = FALSE), detection_rounding)
    dist_m[dist_m > 2] <- NA
    if (length(dist_m) == sum(is.na(dist_m))) {
      out <- data.frame(n.outliers = nrow(outl), n.regular.outliers = 0, 
                        regular.distance = NA)
    }
    else {
      dist_m[dist_m < 10^(-detection_rounding)] <- NA
      dists <- c(dist_m)
      dists <- sort(table(dists))
      dist_m <- as.matrix(dist_m)
      dist_m[row(dist_m) <= col(dist_m)] <- NA
      com_dist <- as.numeric(names(which(dists == max(dists))))
      sel <- which(dist_m %in% com_dist[1])
      sel <- unique(arrayInd(sel, dim(dist_m)))
      reg_outl <- cbind(outl[sel[, 1], 2], outl[sel[, 
                                                    2], 2])
      if (graphs) {
        y0 <- max(x$gamma)
        if (y0 < 0.3 | y0 > 2) {
          y1 <- max(x$gamma) - max(x$gamma)/nrow(reg_outl)
        }
        else {
          y1 <- max(x$gamma) - 0.1
        }
        ys <- seq(y0, y1, by = -((y0 - y1)/(nrow(reg_outl) - 
                                              1)))
        segments(x0 = reg_outl[, 1], x1 = reg_outl[, 
                                                   2], y0 = ys, y1 = ys, col = "red")
      }
      out <- data.frame(n.outliers = nrow(outl), n.regular.outliers = nrow(reg_outl), 
                        regular.distance = com_dist[1])
    }
  }
  return(out)
}

# ras_create Non-exported function from CoordinateCleaner ====
#' @noRd
ras_create <- function(x, lat, lon,  thinning_res){
  # get data extend
  ex <- terra::ext(terra::vect(x[, c(lon, lat)], 
                               geom = c(lon, lat))) + thinning_res * 2
  
  # check for boundary conditions
  if (ex[1] < -180 | ex[2] > 180 | ex[3] < -90 | ex[4] > 90) {
    warning("fixing raster boundaries, assuming lat/lon projection")
    
    if (ex[1] < -180) {
      ex[1] <- -180
    }
    
    if (ex[2] > 180) {
      ex[2] <- 180
    }
    
    if (ex[3] < -90) {
      ex[3] <- -90
    }
    
    if (ex[4] > 90) {
      ex[4] <- 90
    }
  }
  
  # create raster
  ras <- terra::rast(x = ex, resolution = thinning_res)
  
  # set cell ids
  vals <- seq_len(terra::ncell(ras))
  ras <- terra::setValues(ras, vals)
  
  return(ras)
}
#ras dist====
# internal function from coordinate cleaner
# A function to get the distance between raster midpoints and 
#output a data.frame with the distances and the cell IDs as row and column names for cc_outl
#' @noRd
ras_dist <-  function(x, lat, lon, ras, weights) {
  #x = a data.frame of point coordinates, ras = a raster with cell IDs as layer,
  #weight = logical, shall the distance matrix be weightened by the number of
  #points per cell? assign each point to a raster cell
  pts <- terra::extract(x = ras, 
                        y = terra::vect(x[, c(lon, lat)],
                                        geom = c(lon, lat),
                                        crs = ras))
  
  # convert to data.frame
  midp <- data.frame(terra::as.points(ras), 
                     terra::xyFromCell(ras, 1:terra::ncell(ras)))
  
  # retain only cells that contain points
  midp <- midp[midp$lyr.1 %in% unique(pts$lyr.1), , drop = FALSE]
  
  # order
  midp <- midp[match(unique(pts$lyr.1), midp$lyr.1), , drop = FALSE]
  
  
  # calculate geospheric distance between raster cells with points
  dist <- geosphere::distm(midp[, c("x", "y")], 
                           fun = geosphere::distHaversine) / 1000
  
  # set rownames and colnames to cell IDs
  dist <- as.data.frame(dist, row.names = as.integer(midp$lyr.1))
  names(dist) <- midp$lyr.1
  
  if (weights) {
    # approximate within cell distance as half 
    # the cell size, assumin 1 deg = 100km
    # this is crude, but doesn't really matter
    dist[dist == 0] <- 100 * mean(terra::res(ras)) / 2
    
    # weight matrix to account for the number of points per cell
    ## the number of points in each cell
    cou <- table(pts$lyr.1)
    
    ## order
    cou <- cou[match(unique(pts$lyr.1), names(cou))]
    
    # weight matrix, representing the number of distances between or within the cellse (points cell 1 * points cell 2)
    wm <- outer(cou, cou)
    
    # multiply matrix elements to get weightend sum
    dist <- round(dist * wm, 0)
    
    dist <- list(pts = pts, dist = dist, wm = wm)
  } else {
    # set diagonale to NA, so it does not influence the mean
    dist[dist == 0] <- NA
    
    dist <- list(pts = pts, dist = dist)
  }
  
  return(dist)
}

#### BIOGEO USED IMPORTED FUNCTIONS ====
# This is done due to the recent archiving of the pkg

# Bioegeo outliers  ====
#' @title Centroid detection function
#' @description Calculates the outliers using the reverse jackknife procedure (see rjack) and boxplot statistics (using boxplot.stats). 
#' If the value lies 1.5 times beyond the length of the box in the boxplot then it is considered to be an outlier. 
#' @details This is a method and functions first applied in the biogeo package
#' @param rid \emph{numeric}. Row identifier
#' @param species \emph{numeric}. Column of species names
#' @param dups \emph{numeric}. a column of zeros and ones, where ones indicate duplicates
#' @param ev \emph{numeric}. values of the environmental variable
#' @keywords internal
#' @family analysis
#' @author Mark Robertson
#' @importFrom grDevices boxplot.stats
#' @returns a matrix with the test outputs for outlier detection (1/0)
#' @export
#' @examples 
#' rid<-1:20
#' species<-rep("Species A",20)
#' dups=rep(0,20)
#' ev<-c(rnorm(19,mean=20,sd=1),40)
#' a<-biogeo_outliers(rid, species, dups, ev)
biogeo_outliers <- function (rid, species, dups, ev) {
  uspp <- unique(species)
  nr <- length(species)
  ee <- rep(0, nr)
  ee2 <- rep(0, nr)
  for (j in 1:length(uspp)) {
    spp <- uspp[j]
    fsp <- which(species == spp & dups == 0 & !is.na(ev))
    if (length(fsp) >= 10) {
      xc <- ev[fsp]
      ri <- rid[fsp]
      b1 <- boxplot.stats(xc, coef = 1.5)
      xr <- range(b1$stats)
      fe <- which(xc > xr[2] | xc < xr[1])
      ff1 <- ri[fe]
      fe2 <- biogeo_rjack(xc)
      ff2 <- ri[fe2]
      ee[ff1] <- 1
      ee2[ff2] <- 1
    }
  }
  out <- cbind(ee, ee2)
  return(out)
}

# bioegeo rjack ====
#' @noRd
biogeo_rjack <- function (d) {
  xx <- d
  d <- unique(d)
  rng <- diff(range(d))
  mx <- mean(d)
  n <- length(d)
  n1 <- n - 1
  t1 <- (0.95 * sqrt(n)) + 0.2
  x <- sort(d)
  y <- rep(0, n1)
  for (i in 1:n1) {
    x1 <- x[i + 1]
    if (x[i] < mx) {
      y[i] <- (x1 - x[i]) * (mx - x[i])
    }
    else {
      y[i] <- (x1 - x[i]) * (x1 - mx)
    }
  }
  my <- mean(y)
  z <- y/(sqrt(sum((y - my)^2)/n1))
  out <- rep(0, length(xx))
  if (any(z > t1)) {
    f <- which(z > t1)
    v <- x[f]
    if (v < median(x)) {
      xa <- (xx <= v) * 1
      out <- out + xa
    }
    if (v > median(x)) {
      xb <- (xx >= v) * 1
      out <- out + xb
    }
  }
  else {
    out <- out
  }
  f <- which(out == 1)
}

# bioegeo coord2numeric ====
#' @noRd
biogeo_coord2numeric <- function (xn) 
{
  if (is.factor(xn)) {
    xn <- as.character(xn)
  }
  x1 <- as.numeric(xn)
  return(x1)
}
### SEEMS LIKE THEY HAVE NOT BEEN USED
#' # ne_as_sf  Non-exported function from rnaturalearth ====
#' #' @noRd
#' ne_as_sf <- function (x, returnclass = c("sp", "sf")) 
#' {
#'   returnclass <- match.arg(returnclass)
#'   if (returnclass == "sf") 
#'     sf::st_as_sf(x)
#'   else x
#' }
#' 
#' # check_and_return_datatable  Non-exported function from dataPreparation ====
#' #' @noRd
#' check_and_return_datatable <- function (data_set, data_set_name = "data_set") 
#' {
#'   if (!(data.table::is.data.table(data_set) || is.data.frame(data_set) || 
#'         is.matrix(data_set))) {
#'     stop(paste(data_set_name, "should be a data.table, a data.frame or a matrix."))
#'   }
#'   if (nrow(data_set) < 1) {
#'     stop(paste(data_set_name, "should have at least have 1 line."))
#'   }
#'   if (ncol(data_set) < 1) {
#'     stop(paste(data_set_name, "should have at least have 1 column."))
#'   }
#'   if (!data.table::is.data.table(data_set)) {
#'     if (is.data.frame(data_set)) {
#'       data.table::setDT(data_set)
#'     }
#'     else {
#'       data_set <- data.table::as.data.table(data_set)
#'     }
#'   }
#'   if (length(unique(names(data_set))) < length(names(data_set))) {
#'     warning(paste0(data_set_name, ": has column names in double : you should take care of it.", 
#'                    "I changed them to be unique."))
#'     data.table::setDT(data_set, check.names = TRUE)
#'   }
#'   data_set <- data.table::alloc.col(data_set)
#'   return(data_set)
#' }
#' 
#' 
#' # is.verbose  Non-exported function from dataPreparation ====
#' #' @noRd
#' is.verbose <- function (verbose, function_name = "is.verbose") 
#' {
#'   if (!is.logical(verbose)) {
#'     stop(function_name, " verbose should be logical (TRUE or FALSE).")
#'   }
#' }
#' 
#' 
#' # is_ambiguities  Non-exported function from dataPreparation ====
#' #' @noRd
#' is_ambiguities <- function (ambiguities, function_name) 
#' {
#'   if (!is.character(ambiguities) || !ambiguities %in% c("IGNORE", 
#'                                                         "WARN", "SOLVE")) {
#'     stop(paste0(function_name, ": ambiguities should be either IGNORE, WARN or SOLVE."))
#'   }
#' }
#' 
#' 
#' # real_cols  Non-exported function from dataPreparation ====
#' #' @noRd
#' real_cols <- function (data_set, cols, function_name = "real_cols", types = NULL, 
#'                        verbose = TRUE) 
#' {
#'   if (is.null(cols) || length(cols) == 0) {
#'     return(NULL)
#'   }
#'   if (all(cols == "auto")) {
#'     cols <- names(data_set)
#'   }
#'   else {
#'     cols <- reduce_cols_to_existing_ones(data_set = data_set, 
#'                                          cols = cols, function_name = function_name, verbose = verbose)
#'     verbose <- FALSE
#'   }
#'   if (!is.null(types)) {
#'     cols <- reduce_cols_to_correct_type(cols = cols, data_set = data_set, 
#'                                         function_name = function_name, types = types, verbose = verbose)
#'   }
#'   return(cols)
#' }
#' 
#' 
#' # printl  Non-exported function from dataPreparation ====
#' #' @noRd
#' printl <- function (...) 
#' {
#'   args <- list(...)
#'   print(paste(args, collapse = ""))
#' }
#' 
#' 
#' # reduce_cols_to_existing_ones   Non-exported function from dataPreparation ====
#' #' @noRd
#' reduce_cols_to_existing_ones <- function (data_set, cols, function_name, verbose) 
#' {
#'   error_list <- !cols %in% names(data_set)
#'   if (sum(error_list) > 0) {
#'     if (verbose) {
#'       printl(function_name, ": ", paste0(cols[error_list], 
#'                                          collapse = ", "), " aren't columns of the table, i do nothing for those variables")
#'     }
#'     cols <- cols[!error_list]
#'   }
#'   return(cols)
#' }
#' 
#' 
#' 
#' # reduce_cols_to_correct_type Non-exported function from dataPreparation ====
#' #' @noRd
#' reduce_cols_to_correct_type <- function (cols, data_set, function_name, types, verbose) 
#' {
#'   if (all(types == "date")) {
#'     col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
#'                                                               assertthat::is.date)]
#'   }
#'   else if (all(types == "numeric")) {
#'     col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
#'                                                               is.numeric)]
#'   }
#'   else {
#'     col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
#'                                                               class) %in% types]
#'   }
#'   if (sum(col_is_of_wrong_type) > 0) {
#'     if (verbose) {
#'       printl(function_name, ": ", print(cols[col_is_of_wrong_type], 
#'                                         collapse = ", "), " aren't columns of types ", 
#'              paste(types, collapse = " or "), " i do nothing for those variables.")
#'     }
#'     cols <- cols[!col_is_of_wrong_type]
#'   }
#'   return(cols)
#' }
