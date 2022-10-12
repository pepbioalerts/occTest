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
ras_create <- function (x, lat, lon, thinning_res) 
{
  ex <- raster::extent(sp::SpatialPoints(x[, c(lon, lat)])) + 
    thinning_res * 2
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
  ras <- raster::raster(x = ex, resolution = thinning_res)
  vals <- seq_len(raster::ncell(ras))
  ras <- raster::setValues(ras, vals)
  return(ras)
}

# ras_dist Non-exported function from CoordinateCleaner ====
#' @noRd
ras_dist <- function (x, lat, lon, ras, weights) 
{
  pts <- raster::extract(x = ras, y = sp::SpatialPoints(x[, 
                                                          c(lon, lat)]))
  midp <- data.frame(raster::rasterToPoints(ras))
  midp <- midp[midp$layer %in% unique(pts), ]
  midp <- midp[match(unique(pts), midp$layer), ]
  dist <- geosphere::distm(midp[, c("x", "y")], fun = geosphere::distHaversine)/1000
  dist <- as.data.frame(dist, row.names = as.integer(midp$layer))
  names(dist) <- midp$layer
  if (weights) {
    dist[dist == 0] <- 100 * mean(raster::res(ras))/2
    cou <- table(pts)
    cou <- cou[match(unique(pts), names(cou))]
    wm <- outer(cou, cou)
    dist <- round(dist * wm, 0)
    dist <- list(pts = pts, dist = dist, wm = wm)
  }
  else {
    dist[dist == 0] <- NA
    dist <- list(pts = pts, dist = dist)
  }
  return(dist)
}

# .download  Non-exported function from raster ====
#' @noRd
#' setting warn to -1 by default to avoid long download messages from raster package that blur msgs of the main function
.download <- function (aurl, filename, w = -1) 
{
  fn <- paste(tempfile(), ".download", sep = "")
  res <- utils::download.file(url = aurl, destfile = fn, quiet = FALSE, 
                              mode = "wb", cacheOK = TRUE)
  if (res == 0) {
    # w <- getOption("warn")
    # on.exit(options(warn = w))
    # options(warn = -1)
    options(warn= w)
    if (!file.rename(fn, filename)) {
      file.copy(fn, filename)
      file.remove(fn)
    }
  }
  else {
    stop("could not download the file")
  }
}


# ne_as_sf  Non-exported function from rnaturalearth ====
#' @noRd
ne_as_sf <- function (x, returnclass = c("sp", "sf")) 
{
  returnclass <- match.arg(returnclass)
  if (returnclass == "sf") 
    st_as_sf(x)
  else x
}



# check_and_return_datatable  Non-exported function from dataPreparation ====
#' @noRd
check_and_return_datatable <- function (data_set, data_set_name = "data_set") 
{
  if (!(data.table::is.data.table(data_set) || is.data.frame(data_set) || 
        is.matrix(data_set))) {
    stop(paste(data_set_name, "should be a data.table, a data.frame or a matrix."))
  }
  if (nrow(data_set) < 1) {
    stop(paste(data_set_name, "should have at least have 1 line."))
  }
  if (ncol(data_set) < 1) {
    stop(paste(data_set_name, "should have at least have 1 column."))
  }
  if (!data.table::is.data.table(data_set)) {
    if (is.data.frame(data_set)) {
      data.table::setDT(data_set)
    }
    else {
      data_set <- data.table::as.data.table(data_set)
    }
  }
  if (length(unique(names(data_set))) < length(names(data_set))) {
    warning(paste0(data_set_name, ": has column names in double : you should take care of it.", 
                   "I changed them to be unique."))
    data.table::setDT(data_set, check.names = TRUE)
  }
  data_set <- data.table::alloc.col(data_set)
  return(data_set)
}


# is.verbose  Non-exported function from dataPreparation ====
#' @noRd
is.verbose <- function (verbose, function_name = "is.verbose") 
{
  if (!is.logical(verbose)) {
    stop(function_name, " verbose should be logical (TRUE or FALSE).")
  }
}


# is_ambiguities  Non-exported function from dataPreparation ====
#' @noRd
is_ambiguities <- function (ambiguities, function_name) 
{
  if (!is.character(ambiguities) || !ambiguities %in% c("IGNORE", 
                                                        "WARN", "SOLVE")) {
    stop(paste0(function_name, ": ambiguities should be either IGNORE, WARN or SOLVE."))
  }
}


# real_cols  Non-exported function from dataPreparation ====
#' @noRd
real_cols <- function (data_set, cols, function_name = "real_cols", types = NULL, 
                       verbose = TRUE) 
{
  if (is.null(cols) || length(cols) == 0) {
    return(NULL)
  }
  if (all(cols == "auto")) {
    cols <- names(data_set)
  }
  else {
    cols <- reduce_cols_to_existing_ones(data_set = data_set, 
                                         cols = cols, function_name = function_name, verbose = verbose)
    verbose <- FALSE
  }
  if (!is.null(types)) {
    cols <- reduce_cols_to_correct_type(cols = cols, data_set = data_set, 
                                        function_name = function_name, types = types, verbose = verbose)
  }
  return(cols)
}


# printl  Non-exported function from dataPreparation ====
#' @noRd
printl <- function (...) 
{
  args <- list(...)
  print(paste(args, collapse = ""))
}


# reduce_cols_to_existing_ones   Non-exported function from dataPreparation ====
#' @noRd
reduce_cols_to_existing_ones <- function (data_set, cols, function_name, verbose) 
{
  error_list <- !cols %in% names(data_set)
  if (sum(error_list) > 0) {
    if (verbose) {
      printl(function_name, ": ", paste0(cols[error_list], 
                                         collapse = ", "), " aren't columns of the table, i do nothing for those variables")
    }
    cols <- cols[!error_list]
  }
  return(cols)
}



# reduce_cols_to_correct_type Non-exported function from dataPreparation ====
#' @noRd
reduce_cols_to_correct_type <- function (cols, data_set, function_name, types, verbose) 
{
  if (all(types == "date")) {
    col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
                                                              assertthat::is.date)]
  }
  else if (all(types == "numeric")) {
    col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
                                                              is.numeric)]
  }
  else {
    col_is_of_wrong_type <- !cols %in% names(data_set)[sapply(data_set, 
                                                              class) %in% types]
  }
  if (sum(col_is_of_wrong_type) > 0) {
    if (verbose) {
      printl(function_name, ": ", print(cols[col_is_of_wrong_type], 
                                        collapse = ", "), " aren't columns of types ", 
             paste(types, collapse = " or "), " i do nothing for those variables.")
    }
    cols <- cols[!col_is_of_wrong_type]
  }
  return(cols)
}
