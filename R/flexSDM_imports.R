### FUNCTIONS IMPORTED FROM FLEXSDM

#' Integration of outliers detection methods in environmental space using flexsdm pkg
#'
#' @description This function is a  modified function from flexsdm package for better integration to occTest. 
#' From flexsdm. This function performs different methods for detecting outliers in species
#' distribution data based on the environmental conditions of occurrences. Some methods need
#' presence and absence data (e.g. Two-class Support Vector Machine and Random Forest) while other
#' only use presences (e.g. Reverse Jackknife, Box-plot, and Random Forest outliers) .
#' Outlier detection can be a useful procedure in occurrence data cleaning (Chapman 2005, Liu et al., 2018).
#'
#' @param data data.frame or tibble with presence (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data.
#' @param y character. Column name with latitude data.
#' @param id character. Column name with row id. Each row (record) must have its
#' own unique code.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param env_layer SpatRaster. Raster with environmental variables
#'
#' @details
#' From flexsdm: This function will apply outliers detection methods to occurrence data.
#' Box-plot and Reverse Jackknife method will test outliers for each variable individually, if an
#' occurrence behaves as an outlier for at least one variable it will be highlighted as an outlier.
#' If the user uses only presence data, Support Vector Machine and Random Forest Methods will not be
#' performed. Support Vector Machine and Random Forest are performed with default
#' hyper-parameter values. In the case of a species with < 7 occurrences, the function
#' will not perform any methods (i.e. the additional columns will have 0 values); nonetheless, it will return a tibble with the additional columns with 0 and 1.
#' For further information about these methods, see Chapman (2005), Liu et al. (2018), and Velazco
#' et al. (2022).
#' @author Santiago J.E. Velazco
#' @seealso \link[flexsdm]{env_outliers}
#' @return A tibble object with the same database used in 'data' argument and with seven additional columns, where 1 and 0 denote that a presence was detected or not as outliers
#' \itemize{
#'   \item .out_bxpt: outliers detected with Box-plot method
#'   \item .out_jack: outliers detected with Reverse Jackknife method
#'   \item .out_svm: outliers detected with Support Vector Machine method
#'   \item .out_rf: outliers detected with Random Forest method
#'   \item .out_rfout: outliers detected with Random Forest Outliers method
#'   \item .out_lof: outliers detected with Local outlier factor  method
#'   \item .out_sum: frequency of a presences records was detected as outliers
#'   based on the previews methods (values between 0 and 6).
#'   }
#'
#' @references
#' \itemize{
#'   \item Chapman, A. D. (2005). Principles and methods of data cleaning: Primary Species and
#'   Species- Occurrence Data. version 1.0. Report for the Global Biodiversity Information
#'   Facility, Copenhagen. p72.  http://www.gbif.org/document/80528
#'   \item Liu, C., White, M., & Newell, G. (2018). Detecting outliers in species distribution
#'   data. Journal of Biogeography, 45(1), 164 - 176. https://doi.org/10.1111/jbi.13122
#'   \item Velazco, S.J.E.; Bedrij, N.A.; Keller, H.A.; Rojas, J.L.; Ribeiro, B.R.; De Marco, P. (2022)
#'   Quantifying the role of protected areas for safeguarding the uses of biodiversity.
#'   Biological Conservation, xx(xx) xx-xx. https://doi.org/10.1016/j.biocon.2022.109525
#'   }
#'
#' @importFrom dplyr select mutate tibble filter pull starts_with bind_rows
#' @importFrom grDevices boxplot.stats
#' @importFrom kernlab ksvm predict
#' @importFrom randomForest randomForest outlier
#' @importFrom Rlof lof
#' @importFrom stats quantile
#' @importFrom terra extract
#' @examples
#' \dontrun{
#' #Example from flexsdm
#' require(dplyr)
#' require(terra)
#' require(ggplot2)
#'
#' # Environmental variables
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Species occurrences
#' data("spp")
#' spp
#' spp1 <- spp %>% dplyr::filter(species == "sp1")
#'
#' somevar[[1]] %>% plot()
#' points(spp1 %>% filter(pr_ab == 1) %>% select(x, y), col = "blue", pch = 19)
#' points(spp1 %>% filter(pr_ab == 0) %>% select(x, y), col = "red", cex = 0.5)
#'
#' spp1 <- spp1 %>% mutate(idd = 1:nrow(spp1))
#'
#' # Detect outliers
#' outs_1 <- env_outliers(
#'   data = spp1,
#'   pr_ab = "pr_ab",
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar
#' )
#'
#' # How many outliers were detected by different methods?
#' out_pa <- outs_1 %>%
#'   dplyr::select(starts_with("."), -.out_sum) %>%
#'   apply(., 2, function(x) sum(x, na.rm = T))
#' out_pa
#'
#' # How many outliers were detected by the sum of different methods?
#' outs_1 %>%
#'   dplyr::group_by(.out_sum) %>%
#'   dplyr::count()
#'
#' # Let explor where are locate records highlighted as outliers
#' outs_1 %>%
#'   dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point(aes(col = factor(.out_sum))) +
#'   facet_wrap(. ~ factor(.out_sum))
#'
#' # Detect outliers only with presences
#' outs_2 <- env_outliers(
#'   data = spp1 %>% dplyr::filter(pr_ab == 1),
#'   pr_ab = "pr_ab",
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar
#' )
#'
#' # How many outliers were detected by different methods
#' out_p <- outs_2 %>%
#'   dplyr::select(starts_with("."), -.out_sum) %>%
#'   apply(., 2, function(x) sum(x, na.rm = T))
#'
#' # How many outliers were detected by the sum of different methods?
#' outs_2 %>%
#'   dplyr::group_by(.out_sum) %>%
#'   dplyr::count()
#'
#' # Let explor where are locate records highlighted as outliers
#' outs_2 %>%
#'   dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point(aes(col = factor(.out_sum))) +
#'   facet_wrap(. ~ factor(.out_sum))
#'
#'
#' # Comparison of function outputs when using it with
#' # presences-absences or only presences data.
#'
#' bind_rows(out_p, out_pa)
#' # Because the second case only were used presences, outliers methods
#' # based in Random Forest (.out_rf) and Support Vector Machines (.out_svm)
#' # were not performed.
#' }
env_outliers <- function(data, x, y, pr_ab, id, env_layer) {
  . <- NULL
  # Select columns and rename them
  data0 <- data
  data <- data[, c(id, x, y, pr_ab)]
  names(data) <- c("id", "x", "y", "pr_ab")
  
  # Convert data to tibble object
  data <- data %>% tibble()
  var <- names(env_layer)
  
  out_list <- list()
  occ_sp_01 <- data %>%
    dplyr::select(x, y, id, pr_ab)
  
  occ_sp_01 <-
    occ_sp_01 %>% dplyr::mutate(
      .out_bxpt = 0,
      .out_jack = 0,
      .out_svm = 0,
      .out_rf = 0,
      .out_rfout = 0,
      .out_lof = 0,
      .out_sum = 0
    )
  
  sp_env_01 <-
    terra::extract(
      env_layer,
      terra::vect(occ_sp_01[c("x", "y")] %>%
                    dplyr::rename(lon = x, lat = y))
    )[-1] %>%
    data.frame() %>%
    dplyr::tibble(id = occ_sp_01$id, pr_ab = occ_sp_01$pr_ab, .)
  
  # Remove NAs
  complete_vec <- stats::complete.cases(sp_env_01)
  if (sum(!complete_vec) > 0) {
    message(
      sum(!complete_vec),
      " rows were excluded from database because NAs were found"
    )
    occ_sp_01 <- occ_sp_01 %>% dplyr::filter(complete_vec)
    sp_env_01 <- sp_env_01 %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)
  
  sp_env_1 <- sp_env_01 %>% dplyr::filter(pr_ab == 1)
  occ_sp_01 <- occ_sp_01 %>% dplyr::select(-x, -y, -pr_ab)
  
  p01 <- unique(sp_env_01$pr_ab) # vector for testing presence and absence
  
  if (nrow(sp_env_1) > 6) {
    #### Method based on Boxplot and Reverse Jackknife ####
    l_box <- matrix(0, nrow = nrow(sp_env_1), length(var))
    l_jackk <- matrix(0, nrow = nrow(sp_env_1), length(var))
    for (ii in 1:length(var)) {
      xc <- sp_env_1 %>%
        data.frame() %>%
        dplyr::pull(var[ii])
      xr <-
        grDevices::boxplot.stats(xc, coef = 1.5)$stats %>% range()
      fe <- which((xc > xr[2] | xc < xr[1]))
      fe2 <- rev_jack(v = xc) # reverse Reverse Jackknife
      l_box[fe, ii] <- 1
      l_jackk[fe2, ii] <- 1
    }
    
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_bxpt"] <-
      ifelse(rowSums(l_box) > 0, 1, 0)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_jack"] <-
      ifelse(rowSums(l_jackk) > 0, 1, 0)
    
    #### 	Two-class Support Vector Machine (presences absences) ####
    if (all(c(0, 1) %in% p01)) {
      sv <-
        kernlab::ksvm(
          pr_ab ~ .,
          data = sp_env_01[-1] %>% dplyr::mutate(pr_ab = factor(pr_ab)),
          type = "C-bsvc",
          kernel = "rbfdot",
          kpar = list(sigma = 0.1),
          C = 10,
          prob.model = TRUE
        )
      psv2 <-
        kernlab::predict(sv, sp_env_1, type = "probabilities")[, 2] # prediction for presences
      psv2 <- 1 - psv2 # outlierness
      
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_svm"] <-
        as.integer(psv2 > stats::quantile(psv2, probs = seq(0, 1, 0.05))[20], na.rm = TRUE)
      rm(sv)
      
      #### Random Forest ####
      rf <-
        randomForest::randomForest(
          pr_ab ~ .,
          data = data.frame(sp_env_01[-1]) %>% dplyr::mutate(pr_ab = factor(pr_ab)),
          ntree = 2000
        )
      prd2 <- stats::predict(rf, sp_env_1[-1], "prob")[, 2]
      prd2 <- 1 - prd2 # outlierness
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rf"] <-
        as.integer(prd2 > stats::quantile(prd2, probs = seq(0, 1, 0.05))[20], na.rm = TRUE)
      rm(rf)
    }
    
    #### Random forest - outlier ####
    rf <-
      randomForest::randomForest(
        sp_env_1[-1] %>% dplyr::mutate(
          pr_ab =
            factor(pr_ab)
        ),
        ntree = 2000,
        proximity = TRUE
      )
    ot <- randomForest::outlier(rf)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rfout"] <-
      as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    rm(rf)
    
    #### Local outliers factor ####
    if (nrow(sp_env_1) < 15) {
      ot <- rep(NA, nrow(sp_env_1))
      wii = 2
      while (any(is.na(ot))) {
        ot <- Rlof::lof(sp_env_1[-1], k = wii, cores = 1)
        wii <- wii + 1
      }
    } else {
      ot <- rep(NA, nrow(sp_env_1))
      wii=15
      while(any(is.na(ot))){
        ot <- Rlof::lof(sp_env_1[-1], k = ifelse(wii>= nrow(sp_env_1),  nrow(sp_env_1)-1, wii), cores = 1)
        wii <- wii + 5
      }
    }
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_lof"] <-
      as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    rm(ot)
    
    # Summary outliers
    occ_sp_01[, ".out_sum"] <-
      occ_sp_01 %>%
      dplyr::select(dplyr::starts_with(".")) %>%
      rowSums()
    out_list <- occ_sp_01
  } else {
    out_list <- occ_sp_01
  }
  
  out <- dplyr::bind_rows(out_list)
  cn <- "id"
  names(cn) <- id
  out <- dplyr::left_join(data0, out, by = cn)
  return(out)
}


## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##


#' pre_tr_te
#'
#' @noRd
pre_tr_te <- function(data, p_names, h) {
  train <- list()
  test <- list()
  
  if (any(c("train", "train-test", "test")
          %in%
          unique(data[, p_names[h]]))) {
    np2 <- 1
    filt <- grepl("train", data[, p_names[h]])
    train[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])
    filt <- grepl("test", data[, p_names[h]])
    test[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])
  } else {
    np2 <- max(data[p_names[h]])
    
    for (i in 1:np2) {
      train[[i]] <- data[data[p_names[h]] != i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
      
      test[[i]] <- data[data[p_names[h]] == i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
    }
  }
  return(list(train = train, test = test, np2 = np2))
}



# Inverse bioclim
#'
#' @noRd
#'
bio <- function(data, env_layer) {
  . <- NULL
  if (class(data)[1] != "data.frame") {
    data <- data.frame(data)
  }
  if (!methods::is(env_layer, "SpatRaster")) {
    env_layer <- terra::rast(env_layer)
  }
  
  data <- stats::na.omit(data)
  
  result <- env_layer[[1]]
  result[] <- NA
  
  minv <- apply(data, 2, min)
  maxv <- apply(data, 2, max)
  vnames <- names(data)
  
  data_2 <- data %>%
    stats::na.omit() %>%
    apply(., 2, sort) %>%
    data.frame()
  
  rnk <- function(x, y) {
    b <- apply(y, 1, FUN = function(z) sum(x < z))
    t <- apply(y, 1, FUN = function(z) sum(x == z))
    r <- (b + 0.5 * t) / length(x)
    i <- which(r > 0.5)
    r[i] <- 1 - r[i]
    r * 2
  }
  
  var_df <- terra::as.data.frame(env_layer)
  var_df <- stats::na.omit(var_df)
  
  k <- (apply(t(var_df) >= minv, 2, all) &
          apply(t(var_df) <= maxv, 2, all))
  
  for (j in vnames) {
    var_df[k, j] <- rnk(
      data_2[, j],
      var_df[k, j, drop = FALSE]
    )
  }
  var_df[!k, ] <- 0
  res <- apply(var_df, 1, min)
  result[as.numeric(names(res))] <- res
  return(result)
}

inv_bio <- function(e, p) {
  if (!methods::is(e, "SpatRaster")) {
    e <- terra::rast(e)
  }
  r <- bio(data = terra::extract(e, p)[-1], env_layer = e)
  r <- (r - terra::minmax(r)[1]) /
    (terra::minmax(r)[2] - terra::minmax(r)[1])
  r <- r <= 0.01 # environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}


#' Inverse geo
#'
#' @noRd
#'
inv_geo <- function(e, p, d) {
  colnames(p) <- c("x", "y")
  p <- terra::vect(p, geom = c("x", "y"), crs = terra::crs(e))
  b <- terra::buffer(p, width = d)
  b <- terra::rasterize(b, e, background = 0)
  e <- terra::mask(e, b, maskvalues = 1)
  return(e)
}

#' Boyce
#'
#' @description This function calculate Boyce index performance metric. Codes were adapted from
#' enmSdm package.
#'
#' @noRd
boyce <- function(pres,
                  contrast,
                  n_bins = 101,
                  n_width = 0.1) {
  lowest <- min(c(pres, contrast), na.rm = TRUE)
  highest <- max(c(pres, contrast), na.rm = TRUE) + .Machine$double.eps
  window_width <- n_width * (highest - lowest)
  
  lows <- seq(lowest, highest - window_width, length.out = n_bins)
  highs <- seq(lowest + window_width + .Machine$double.eps, highest, length.out = n_bins)
  
  ## initiate variables to store predicted/expected (P/E) values
  freq_pres <- NA
  freq_contrast <- NA
  
  # tally proportion of test presences/background in each class
  for (i in 1:n_bins) {
    # number of presence predictions in a class
    freq_pres[i] <-
      sum(pres >= lows[i] & pres < highs[i], na.rm = TRUE)
    
    # number of background predictions in this class
    freq_contrast[i] <-
      sum(contrast >= lows[i] & contrast < highs[i], na.rm = TRUE)
  }
  
  # mean bin prediction
  mean_pred <- rowMeans(cbind(lows, highs))
  
  # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
  if (any(freq_pres > 0 & freq_contrast == 0)) {
    small_value <- 0.5
    freq_contrast[freq_pres > 0 & freq_contrast == 0] <- small_value
  }
  
  # remove classes with 0 presence frequency
  if (any(freq_pres == 0)) {
    zeros <- which(freq_pres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }
  
  # remove classes with 0 background frequency
  if (any(0 %in% freq_contrast)) {
    zeros <- which(freq_pres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }
  
  P <- freq_pres / length(pres)
  E <- freq_contrast / length(contrast)
  PE <- P / E
  
  # remove NAs
  rm_nas <- stats::complete.cases(data.frame(mean_pred, PE))
  # mean_pred <- mean_pred[rm_nas]
  # PE <- PE[rm_nas]
  
  # calculate Boyce index
  result <- stats::cor(
    x = ifelse(is.na(mean_pred), 0, mean_pred),
    y = ifelse(is.na(PE), 0, PE), method = "spearman"
  )
  return(result)
}

# maxnet:::predict.maxnet()
#' Predict maxnet
#' @importFrom stats model.matrix
#' @importFrom stats formula
#' @noRd
predict_maxnet <- function(object, newdata, clamp = TRUE, type = c("link", "exponential", "cloglog", "logistic"), ...) {
  categoricalval <- function(x, category) {
    ifelse(x == category, 1, 0)
  }
  thresholdval <- function (x, knot)
  {
    ifelse(x >= knot, 1, 0)
  }
  hingeval <- function(x, min, max) {
    pmin(1, pmax(0, (x - min) / (max - min)))
  }
  
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }
  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  mm <- model.matrix(f, data.frame(newdata))
  if (clamp)
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  link <- (mm %*% object$betas) + object$alpha
  type <- match.arg(type)
  
  if (type == "link"){
    return(link)
  }
  if (type == "exponential"){
    return(exp(link))
  }
  if (type == "cloglog"){
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic"){
    return(1/(1 + exp(-object$entropy - link)))
  }
}

#' Outliers with Reverse Jackknife
#' @importFrom stats median
#' @noRd
#'
rev_jack <- function(v) {
  v2 <- v
  v <- unique(v)
  lgh <- length(v) - 1
  t1 <- (0.95 * sqrt(length(v))) + 0.2
  x <- sort(v)
  y <- rep(0, lgh)
  for (i in seq_len(lgh)) {
    x1 <- x[i + 1]
    if (x[i] < mean(v)) {
      y[i] <- (x1 - x[i]) * (mean(v) - x[i])
    } else {
      y[i] <- (x1 - x[i]) * (x1 - mean(v))
    }
  }
  my <- mean(y)
  z <- y / (sqrt(sum((y - my)^2) / lgh))
  out <- rep(0, length(v2))
  if (any(z > t1)) {
    f <- which(z > t1)
    v <- x[f]
    for (vt in v) {
      if (vt < median(x)) {
        xa <- (v2 <= vt) * 1
        out <- out + xa
      }
      if (vt > median(x)) {
        xb <- (v2 >= vt) * 1
        out <- out + xb
      }
    }
  } else {
    out <- out
  }
  #old flexsdm code that breaks
  # if (any(z > t1)) {
  #   f <- which(z > t1)
  #   v <- x[f]
  #   if (v < median(x)) {
  #     xa <- (v2 <= v) * 1
  #     out <- out + xa
  #   }
  #   if (v > median(x)) {
  #     xb <- (v2 >= v) * 1
  #     out <- out + xb
  #   }
  # } else {
  #   out <- out
  # }
  return(which(out == 1))
}

#' Calculate amount of data for each training dataset in a given partition
#'
#' @noRd
#'
n_training <- function(data, partition) {
  . <- partt <- NULL
  if (any(c("train", "train-test", "test")
          %in%
          (data %>%
           dplyr::select(dplyr::starts_with({{partition}})) %>%
           dplyr::pull() %>%
           unique()))) {
    nn_part <- data %>%
      dplyr::select(dplyr::starts_with({{partition}})) %>%
      apply(., 2, table) %>%
      data.frame()
    nn_part <- nn_part %>% dplyr::mutate(partt = rownames(nn_part))
    nn_part$partt[grepl("train", nn_part$partt)] <- "train"
    nn_part <- nn_part %>%
      dplyr::filter(partt == "train") %>%
      dplyr::select(-partt)
    nn_part <- colSums(nn_part)
  } else {
    data <- data %>%
      dplyr::select(dplyr::starts_with({{partition}}))
    
    nn_part <- list()
    for(ppp in 1:ncol(data)){
      nn_part[[ppp]] <- data %>% dplyr::pull(ppp) %>% table() %>% c()
      sm <- nn_part[[ppp]] %>% sum()
      nn_part[[ppp]] <- sapply(nn_part[[ppp]], function(x) sum(sm-x))
    }
    nn_part <- unlist(nn_part)
  }
  return(nn_part)
}

#' Calculate number of coefficient for gam models
#'
#' @noRd
#'
n_coefficients <- function(data, predictors, predictors_f = NULL, k = 10){
  data <- data.frame(data)
  if(k<0){
    k=10
  }
  if(!is.null(predictors_f)){
    n_levels <- rep(NA, length(predictors_f))
    for(fff in 1:length(predictors_f)){
      n_levels[fff] <- unique(data[, predictors_f]) %>% stats::na.omit() %>% length()
    }
    n_levels <- sum(n_levels)
  }else{
    n_levels <- 0
  }
  n <- (k - 1) * length(predictors) + n_levels
  return(n)
}

#' Euclidean distance for extrapolation
#'
#' @noRd
#'
euc_dist <- function(x, y) {
  if (!methods::is(x, "matrix")) {
    x <- as.matrix(x)
  }
  if (!methods::is(y, "matrix")) {
    y <- as.matrix(y)
  }
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (ii in 1:nrow(y)) {
    result[, ii] <- sqrt(colSums((t(x) - y[ii, ]) ^ 2))
  }
  rownames(result) <- rownames(x)
  colnames(result) <- rownames(y)
  return(result)
}

#' Moran I, based on ape package
#' @importFrom stats sd
#' @noRd
#'
morani <- function(x, weight, na.rm = FALSE, scaled = TRUE) {
  if (dim(weight)[1] != dim(weight)[2]) {
    stop("'weight' must be a square matrix")
  }
  n <- length(x)
  if (dim(weight)[1] != n) {
    stop("'weight' must have as many rows as observations in 'x'")
  }
  ei <- -1 / (n - 1)
  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      n <- length(x)
      weight <- weight[!nas, !nas]
    } else {
      warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
      return(NA)
    }
  }
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight / rs
  s <- sum(weight)
  m <- mean(x)
  y <- x - m
  cv <- sum(weight * y %o% y)
  v <- sum(y^2)
  res <- (n / s) * (cv / v)
  if (scaled) {
    imax <- (n / s) * (sd(rowSums(weight) * y) / sqrt(v / (n - 1)))
    res <- res / imax
  }
  
  return(res)
}

#' Mahalanobis distance
#' @importFrom stats mahalanobis
#' @noRd
#'
mah_dist <- function(x, y, cov) {
  if (!methods::is(x, "matrix")) {
    x <- as.matrix(x)
  }
  if (!methods::is(y, "matrix")) {
    y <- as.matrix(y)
  }
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (ii in 1:nrow(y)) {
    # root square of squared Mahalanobis distance
    result[, ii] <- sqrt(mahalanobis(x = x, center = y[ii, ], cov = cov))
  }
  rownames(result) <- rownames(x)
  colnames(result) <- rownames(y)
  return(result)
}
