#### FUNCTIONS FOR ANALYZING DATA IN THE WORKFLOW

# filterMissing ====
#' @title Check for missing coordinates
#' @description checks for missing coordinates in the occurrence dataframe
#' @param df dat.frame of species occurrences
#' @param xf character. The field in the data.frame containing the x coordinates
#' @param yf character. The field in the data.frame containing the y coordinates
#' @param verbose logical. Print messages? Default FALSE
#' @return List with two dataframes: stay = coordinates missing and continue = occurrence that you can retain for further analysis
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples 
#' ds <- data.frame (x=c(runif (n = 100),NA),y=c(runif (n = 100),1000))
#' filterMissing(ds,xf='x',yf='y')
#' @export
filterMissing <- function (df, xf , yf , verbose=FALSE){
 coordIssues_coordMissing_value = !stats::complete.cases  (df[,c(xf,yf)])
 df$coordIssues_coordMissing_value = coordIssues_coordMissing_value
 df$coordIssues_coordMissing_test = as.logical (coordIssues_coordMissing_value)
 df.out <- df [coordIssues_coordMissing_value,]
 df.continue <- df [!coordIssues_coordMissing_value,]
 return (list (stay=df.out,continue=df.continue))
}

# duplicatesexcludeAnalysis ====
#' @title Duplicated records
#' @description checks for duplicated coordinates in the occurrence data frame
#' @details it differentiates the exact duplicates and the duplicates for a occurrences falling in the same pixel
#' @param df Data.frame of species occurrences
#' @param xf the field in the data frame containing the x coordinators
#' @param yf the field in the data frame containing the y coordinates
#' @param resolution.in.minutes the resolution of environmental data used, specified in minutes
#' @param raster.grid An optional raster grid
#' @param verbose logical. Print messages? Default FALSE
#' @return list of three components: Dups.Exact = exact duplicate records, Dups.Grid= duplicates within environmental gridcell, continue = dataframe with good records
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @family analysis
#' @examples 
#' ds <- data.frame (x=c(runif (n = 100),1000),y=c(runif (n = 100),1000),Reason=NA)
#' duplicatesexcludeAnalysis(ds,xf='x',yf='y',resolution.in.minutes=60)
#' @export
duplicatesexcludeAnalysis <- function (df, xf, yf,
                                        resolution.in.minutes,
                                        raster.grid=NULL, verbose=FALSE){
  #exact duplicates
  duplicates_dupExact_value <- duplicated (x = df[,c(xf,yf)]) *1
  duplicates_dupExact_test <- as.logical (duplicates_dupExact_value)
  exact.dups = duplicates_dupExact_value
  df$Exclude <- exact.dups
  df$Reason [which (exact.dups==1)] <-'Exact Duplicated coordinates'
  df$duplicates_dupExact_value = duplicates_dupExact_value
  df$duplicates_dupExact_test = duplicates_dupExact_test
  df.exact.dups <- df[which(df$Exclude==1),]
  #update df
  df <- df[which(df$Exclude==0),]
  #build the grid for the duplicates
  if(!is.null(raster.grid)){
    rst <- raster.grid
    }
  if(is.null(raster.grid)) {
    rst <- terra::rast(extent=c(-180, 180, -90,  90),
                  res = resolution.in.minutes/60)
  }
  #get dups in gridcell
  spp <- as.character(unique(df$Species))
  xy <- data.frame(biogeo::coord2numeric(df[,xf]),
                   biogeo::coord2numeric(df[,yf]))
  cid <- terra::cellFromXY(rst, xy)
  dups <- (duplicated(cid)) * 1
  df$duplicates_dupGrid_value = dups
  df$duplicates_dupGrid_test = as.logical (dups)
  f1 <- which (dups==1)
  df$Exclude <- dups
  df$Reason[f1] <- "Duplicated--GridCell"

  #indicate duplicates by grid cell (these needs updateing for efficiency)
  df.grid.dups <- df[which(df$Exclude==1),]

  #update df
  df <- df[which(df$Exclude==0),]

  #return
  return (list (Dups.Exact=df.exact.dups, Dups.Grid= df.grid.dups, continue = df))
}
# SeaLand reassignment ====
#' @title SeaLand reassignment
#' @description Reassign coastal coordinates as needed to be in 
#' @details (function inspired in nearestcell in biogeo but modified 
#' @param dat Data.frame of species occurrences 
#' @param rst raster class object.
#' @param xf the field in the data frame containing the x coordinators
#' @param yf the field in the data frame containing the y coordinators
#' @param verbose logical. Print messages? Default FALSE
#' @return list
#' @keywords internal
#' @author Mark Robertson and Veron Visser (original \link[biogeo]{biogeo-package}), Josep M Serra-Diaz (modifs)
#' @family analysis
#' @seealso original code at \link[biogeo]{nearestcell}
.nearestcell3 <- function (dat,
                          rst,
                          xf,
                          yf,
                          verbose=FALSE) {
  dat$irow <- 1:nrow (dat)
  dat$Correction <- as.character(dat$Correction)
  fx <- which(dat$Exclude == 0)
  x1 <- biogeo::coord2numeric(dat[,xf][fx])
  y1 <- biogeo::coord2numeric(dat[,yf][fx])
  datid <- dat$irow[fx]
  dd <- data.frame(x1, y1)
  ce0 <- terra::cellFromXY(rst, dd)
  vals <- terra::extract(rst, dd)
  #if the input raster is a dataframe
  if ('matrix' %in% class(vals)){
    vals <- apply(vals,1,FUN = sum)
    f <- which(is.na(vals))
  } else {
    f <- which(is.na(vals))
  }
  if (length(f) == 0){ if(verbose) print ("There are no missing values") ; return (dat)}
  if (length (f) != 0){
    id  <- datid[f]
    ce1 <- terra::cellFromXY(rst, dd[f, ])
    if (any (is.na(ce1))) {
      if (length(ce1) == 1) {if(verbose) print("Coordinates in sea are out of environmental extent"); return(dat)}
      if(verbose) print("Some coordinates out of raster extent")
    }
    ff  <- which(!is.na(ce1))
    ce3 <- ce1[ff]
    dd2 <- dd[ff, ]
    #to eliminiate
    id2 <- id[ff]
    bb  <- {}
    for (i in 1:length(ce3)) {
      a <- terra::adjacent(rst, ce3[i], directions = 8, pairs = FALSE,
                     include = FALSE)
      a <- as.numeric(a)
      a <- a[!is.na(a)]
      idx <- id2[i]
      # idx <- ce3[i] 
      b <- data.frame(i, a, id2 = idx)
      bb <- rbind(bb, b)
    }
    xy <- terra::xyFromCell(rst, bb$a)
    vals <- terra::extract(rst[[1]], xy) #tanking only the first layer to inform vals
    g <- data.frame(bb, xy, vals)
    if (any(!is.na(vals))) {
      fv <- which(!is.na(vals))
      id <- bb$i[fv]
      uid <- unique(id)
      g1 <- stats::na.omit(g)
      near <- {
      }

      for (j in 1:length(uid)) {
        uj <- uid[j]
        fx <- which(g1$i == uj)
        nce <- g1[fx, ]
        nce2 <- as.matrix(cbind(nce$x, nce$y))
        pts <- dd2[uj, ]
        k1 = terra::vect(matrix (c(pts$x1, pts$y1),ncol=2),
                         crs="+proj=longlat +datum=WGS84")
        k2 = terra::vect(nce2,
                         crs="+proj=longlat +datum=WGS84")
        dst <- terra::distance(k1, k2)
        fm <- which.min(dst)
        nr <- nce[fm, ]
        near <- rbind(near, nr)
      }
      mod <- format(Sys.time(), "%d-%m-%Y %H:%M:%S")
      dat$Modified <- as.character(dat$Modified)
      for (i in 1:nrow(near)) {
        f <- which(dat$irow == near$id2[i])
        dat$x_original[f] <- dat[,xf][f]
        dat$y_original[f] <- dat[,yf][f]
        dat[,xf][f] <- near$x[i]
        dat[,yf][f] <- near$y[i]
        dc <- dat$Correction[f]
        dc <- stringr::str_replace(dc, "[.]", "7")
        dat$Correction[f] <- dc
        dat$Modified[f] <- mod
      }
      moved <- data.frame(ID = near$id2, x = near$x, y = near$y)
      names(moved) <- c('ID',xf,yf)
      #remove irow column
      columnOut= which (names(dat)=='irow')* (-1)
      dat <- dat[,columnOut]
      datx <- list(dat = dat, moved = moved)
      return(datx)
    } else {
      warning("There are no records close enough to the nearest land/sea cells") 
      return(dat)
    }
  }
}
# countryStatusRange analysis ====
#' @title Range analysis
#' @description Identify and filter species records based on countries where species is considered native or alien.
#' @details It returns a list with two elements: stay (records that are excluded/filtered), continue (records that are not excluded from the filtering process)
#' @param df data.frame of species occurrences
#' @param xf character. The field in the data.frame containing the x coordinates
#' @param yf character. The field in the data.frame containing the y coordinates
#' @param .ntv.ctry character. Native country in ISO3 coding
#' @param .inv.ctry character. Invasive country in ISO3 coding
#' @param .c.field character. Country field in the species data.frame (df)
#' @param .points.crs crs for the occurrence data
#' @param .countries.shapefile spatialPolygonDataFrame of political divisions
#' @param cfsf character. Column name of the .country spatialPolygonDataFrame indicaing the country in ISO3 codeing.
#' @param excludeUnknownRanges logical. Should records be filtered if located in countries outside .ntv.ctry or .inv.ctry? Defalut FALSE
#' @param excludeNotmatchCountry Should records be if the reported country is different than the locatoin country? . Default FALSE
#' @param doRangeAnalysis logical. Should range analysis be performed?
#' @param verbose logical. Print messages? Default FALSE
#' @return list
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @family analysis
#' @export
countryStatusRangeAnalysis=function(df,
                                    xf,
                                    yf,
                                    .ntv.ctry,
                                    .inv.ctry,
                                    .c.field,
                                    .points.crs,
                                    .countries.shapefile,
                                    cfsf,
                                    excludeUnknownRanges = FALSE,
                                    excludeNotmatchCountry = FALSE,
                                    doRangeAnalysis=TRUE,
                                    verbose=FALSE) {
  
  

  if (!doRangeAnalysis) {
    df$countryStatusRange_countryCoordinates   <- NA
    df$countryStatusRange_wrongNTV_value   <- NA
    df$countryStatusRange_wrongNTV_test   <- NA
    df$countryStatusRange_wrongCTRY_value  <- NA
    df$countryStatusRange_wrongCTRY_test  <- NA
    df$countryStatusRange_wrongINV_value <- NA
    df$countryStatusRange_wrongINV_test  <- NA
    df$countryStatusRange_score <-NA
    stay = df[-1*1:nrow(df),]
    out <- list (stay=stay,continue = df)
    return (out)
  }
  #initial user information printing
  if (is.null (.ntv.ctry)  ) {if(verbose) print ('INFO: Species with no associated country. We assume all locations are native range')}
  if (is.null (.inv.ctry)  ) {if(verbose) print ('INFO: No invasive country provided. Analysis of invasive ranges not performed')}
  if (is.null (.c.field)   ) {if(verbose) print ('INFO: No info on table of country of registration. Analysis of country recorded vs. coordinates not performed')}
  xydat <- df[,c(xf,yf)]
  #check country of coordinates
    #extract countries of the points with a given spdf object
    if (!is.null(.countries.shapefile)){
      if (is.null(.points.crs)) {
        .points.crs <- sf::st_crs(.countries.shapefile)
      if(verbose) print ('ASSUMING points in projection')
      }
      sp.xydat <- sf::st_as_sf(xydat,coords=c(xf,yf),crs = .points.crs)
      if (sf::st_crs(sp.xydat) != sf::st_crs(.countries.shapefile))
        .countries.shapefile = sf::st_transform(.countries.shapefile,st_crs(sp.xydat))
      sf::sf_use_s2(F)
      overlay.sp.xydat <- suppressMessages (as.data.frame (sf::st_join(sp.xydat,.countries.shapefile)))
      fieldname  <- match(cfsf, names(overlay.sp.xydat))
      country_ext <- overlay.sp.xydat[, fieldname]
    }
    #extract countries of the points with when a country spdf is not provided
    if (is.null(.countries.shapefile)){
      country_ext <-  .coords2country (xydat)
    }
  #inform for whether reported and country in coordinates (extracted country) differ 
  if (is.null (.c.field)) {
    if(verbose) print ("No country reported in occurrence database")
    wrong.ctry.reported <- rep (NA,length(country_ext))

  }
  if (!is.null (.c.field)) {

    #check if ISO3 format in country reported [FUTURE: use GNRS to automatically correct]
    ctry.reported = as.character (df[,.c.field])
    if (all (nchar (ctry.reported)!=3,na.rm=TRUE) | all (is.na(ctry.reported)) ) {
      warning ("the countries in your database do not match ISO3 codes, test of reported vs. match dropped")
      wrong.ctry.reported <- rep (NA,length(country_ext))
    } else {
      #match reported country with extracted country
      wrong.ctry.reported <- (! as.character (df[,.c.field]) %in% as.character(country_ext))  * 1
      nonMissingVals = ifelse (test = is.na (as.character (df[,.c.field])) | as.character (df[,.c.field]) == '',yes = NA,no = TRUE)
      wrong.ctry.reported <- wrong.ctry.reported * nonMissingVals
    }

  }
  #inform for which locations are not in the native range
  if(is.null (.ntv.ctry)) {
    if(verbose) print ('No info of native country. Assuming all records are native')
    wrong.ntv.ctry.xy <- rep (NA,length(country_ext))
  }
  if(!is.null(.ntv.ctry) ) {
    wrong.ntv.ctry.xy <- (! country_ext %in% .ntv.ctry) * 1
    }

  #inform for which locations are not in the invasive range
  if (is.null (.inv.ctry))   {
    if(verbose) print ('No info of invasive countries')
    wrong.inv.ctry.xy <- rep (NA,length(country_ext))
  }
  if (!is.null (.inv.ctry))  {wrong.inv.ctry.xy <- (! country_ext %in% .inv.ctry) * 1}
  #add to df
  df$countryStatusRange_wrongNTV_value   <- wrong.ntv.ctry.xy
  df$countryStatusRange_wrongNTV_test   <- as.logical (wrong.ntv.ctry.xy)
  df$countryStatusRange_wrongCTRY_value  <- wrong.ctry.reported
  df$countryStatusRange_wrongCTRY_test  <- as.logical (wrong.ctry.reported)
  df$countryStatusRange_wrongINV_value   <- wrong.inv.ctry.xy
  df$countryStatusRange_wrongINV_test   <- as.logical (wrong.inv.ctry.xy)
  df$countryStatusRange_countryCoordinates <- country_ext
  
  #Exclusion from subsequent analysis for wrong native and invasive ranges
  #  (depending how much data has been provided)
  if (excludeUnknownRanges == TRUE & !is.null(.ntv.ctry)  & !is.null(.inv.ctry) ){
    total.data <- (!is.na (wrong.ntv.ctry.xy) *1) + (!is.na(wrong.inv.ctry.xy)*1)
    sum.wrong.extracted.ranges <-  mapply (function (x,y) (sum (x,y,na.rm = TRUE)),wrong.ntv.ctry.xy, wrong.inv.ctry.xy)
    exclusion.based.on.wrong <- sum.wrong.extracted.ranges / total.data
    df$Exclude  <- (exclusion.based.on.wrong == 1)*1
    exclude.rows <- which (df$Exclude == 1)
    #add the reason of exclusion
    if (length(exclude.rows) > 0) {df$Reason [exclude.rows] <-"XY in not in invasive or native ranges"}
  }
  if (excludeUnknownRanges == FALSE) {df$Exclude  <- 0}

  # Exclusion based on wrong country reported
  if  (excludeNotmatchCountry == TRUE & !is.null (.c.field)) {
    exclusion.based.on.wrong <- df$Exclude  + (10*wrong.ctry.reported)
    only.wrong.reported.range <- which (exclusion.based.on.wrong == 10 )
    wrong.reported.range.and.informedrange <- which (exclusion.based.on.wrong == 11 )
    df$Reason [only.wrong.reported.range]  <- "Incongruence reported country and coordinates"
    df$Exclude [only.wrong.reported.range] <- 1
    df$Reason [wrong.reported.range.and.informedrange] <-"XY in not in invasive or native ranges;Incongruence reported country and coordinates"
  }

  #output
  df$countryStatusRange_score =  .gimme.score (df %>% dplyr::select(dplyr::starts_with('countryStatusRange')))
  out <- list (stay = df[which(df$Exclude==1),], continue = df[which(df$Exclude!=1),])
  return (out)
}

# Centroid detection  ====
#' @title Centroid detection function
#' @description Identify occurrence records located near centroids. 
#' @details This is a wrapper function that encompasses different methods
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param .c.field character. Country field in the species data.frame (df)
#' @param .ntv.ctry character. ISO3 country codes where species are considered native
#' @param idf character. Column with the taxon observation ID.
#' @param .inv.ctry character. ISO3 country codes where species are considered alien
#' @param .points.crs crs argument. Coordinate reference system
#' @param .r.env raster. Raster of environmental variables considered in the analysis
#' @param .countries.shapefile spatialPolygonDataFrame of political divisions
#' @param cfsf character. Column name of the .country spatialPolygonDataFrame indicating the country in ISO3 coding.
#' @param method character. Vector with the methods to detect centroids
#' @param do logical. Should range analysis be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @return data.frame
#' @keywords internal
#' @family analysis
#' @details Current methods implemented for centroid detection are 'BIEN' (uses iterative procedure with threshold selection of distance to centroid)\cr
#' , and method 'CoordinateCleaner' implementing methods used in \link[CoordinateCleaner]{CoordinateCleaner-package} package. 
#' @author Brian Maitner, Josep M Serra-Diaz, Cory Merow. Functions from CoordinateCleaner implemented by A Zizka (CoordinateCleaner)
#' @seealso \link[CoordinateCleaner]{cc_cap} \link[CoordinateCleaner]{cc_cen} 
#' @export

centroidDetection <- function (df,
                                xf,
                                yf,
                                cf,
                                idf,
                                .ntv.ctry,
                                .inv.ctry,
                                .points.crs,
                                .r.env,
                                .countries.shapefile,
                                cfsf,
                                 method='all',
                                 do = TRUE, verbose=FALSE){
  #table of outputs
  out = data.frame (centroidDetection_BIEN_value=NA,
                    centroidDetection_BIEN_test=NA,
                    centroidDetection_BIEN_comments =NA,
                    centroidDetection_CoordinateCleaner_value=NA,
                    centroidDetection_CoordinateCleaner_test=NA,
                    centroidDetection_CoordinateCleaner_comment=NA,
                    centroidDetection_score=NA
                    )[1:nrow (df),]

  row.names(out) <- NULL

  if(!do) { return (out)}
  #Method BIEN
  if (any(method %in% c('BIEN','all'))){
    dfxy <- df[,c(xf,yf)]
    if (!is.null(cf)) dfxy$country  <- df[,cf]
    if (is.null(cf))  dfxy$country  <- NA
    if ('countryStatusRange_countryCoordinates' %in% names(df))  dfxy$countryCoordinates  <- df$countryStatusRange_countryCoordinates
    occ_sf <- sf::st_as_sf(dfxy,coords=c(xf,yf),crs = st_crs(4326)) 
    #load centroid data
    dest_url = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/long_centroids.combined.rds'
    outFile = paste0(tempdir(),'/centroids.rds')
    if (!file.exists(outFile)) utils::download.file(url=dest_url,destfile = outFile)
    centroid_data =  readRDS (outFile)
    centroid_data$countryISO3 <- countrycode::countrycode(sourcevar = centroid_data$country,origin = 'country.name' ,destination = 'iso3c')
    bien_centroids <- try (centroid_assessment_BIEN (occ_sf,centroid_data),silent=T)
    if (!class(bien_centroids) %in% 'try-error') {
      out$centroidDetection_BIEN_value <- bien_centroids$is_centroid
      out$centroidDetection_BIEN_test <- as.logical (bien_centroids$is_centroid)
      out$centroidDetection_BIEN_comments <- paste('centroid_dist =',bien_centroids$relative_distance)
    }
    
    if (class(bien_centroids) %in% 'try-error') {
      out$centroidDetection_BIEN_value <- NA
      out$centroidDetection_BIEN_test <- NA
      out$centroidDetection_BIEN_comments <- 'error in bien_centroids'
    }
      
  }
  #Method CoordinateCleaner
  if (any(method %in% c('CoordinateCleaner','all'))){

    cc_cap_test = CoordinateCleaner::cc_cap(x = df,lon =xf,lat=yf,value='flagged', verbose = FALSE)
    cc_cap_test = (!cc_cap_test) *1
    cc_cen_test <- CoordinateCleaner::cc_cen (x = df, lon =xf,lat=yf ,value='flagged', verbose = FALSE)
    cc_cen_test <- (!cc_cen_test) * 1
    out$centroidDetection_CoordinateCleaner_value =  cc_cen_test+cc_cen_test
    out$centroidDetection_CoordinateCleaner_test = ifelse ((cc_cap_test==1 | cc_cen_test==1),TRUE,FALSE)
    out$centroidDetection_CoordinateCleaner_comment = ''
    out$centroidDetection_CoordinateCleaner_comment [cc_cap_test==1] = 'capital centroid'
    out$centroidDetection_CoordinateCleaner_comment [cc_cen_test==1] = 'ctry/prov centroid'

  }
  out$centroidDetection_score =  .gimme.score (x = out)
  return (out)

}

# human detection ====
#' @title Hyper Human Influence
#' @description Detect occurrences in heavily human-impacted environments
#' @details It uses several methods to detect records in high human influence records.\cr
#' Current implemented methods are: \cr
#' 'hii' using a raster and a threhold of human influence
#' 'urban' using a layer of urban areas. Method implemented in CoordinateCleaner package. 
#' @param df Data.frame of species occurrences
#' @param xf the field in the data frame containing the x coordinates
#' @param yf the field in the data frame containing the y coordinates
#' @param ras.hii Raster of human influence index
#' @param .th.human.influence threshold of human influence index
#' @param method character. Indicate which tests to use. Default 'all'. See Details
#' @param .points.crs character. Points coordinate projection.
#' @param ras.hii raster. Raster map of human influence index use
#' @param .th.human.influence numeric. Threhold to identify places of high human influence
#' @param .points.crs crs argument for df
#' @param do logical. Should range analysis be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @param output.dir character. Output directory 
#' @return data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr). A Zizka (CoordinateCleaner functions)
#' @seealso \link[CoordinateCleaner]{cc_urb} for CoordinateCleaner functions
#' @family analysis
#' @export
#' @import sf

humanDetection <- function (df,
                            xf,
                            yf,
                            method='all',
                            ras.hii,
                            .th.human.influence,
                            .points.crs,
                            do=TRUE, verbose=FALSE,output.dir){

  out <- data.frame (humanDetection_HumanInfluence_value=NA,
                     humanDetection_HumanInfluence_test=NA,
                     humanDetection_HumanInfluence_comments= NA,
                     humanDetection_UrbanAreas_value=NA,
                     humanDetection_UrbanAreas_test=NA,
                     humanDetection_UrbanAreas_comments=NA,
                     humanDetection_score=NA
                     )[1:nrow (df),]
  row.names(out) <- NULL

  if (!do) {return (out)}

  #start human influence analysis
  xydat <- df[,c(xf,yf)]

  #get species presence
  data.sp.pres <- as.data.frame (xydat)
  data.sp.pres <- sf::st_as_sf(data.sp.pres,coords=c(xf,yf),crs = .points.crs)

  #start human influence index test
  if (any (method %in% c('hii','all'))) {

    class(ras.hii)=='RasterLayer'

    #accomodate projections
    if (is.null (.points.crs)) {sf::st_set_crs(data.sp.pres, value=terra::crs(ras.hii)) 
      warning ('ASSUMING Points and HumanRaster data have the same projection')
      }
    if (!is.null (.points.crs)) {sf::st_set_crs(data.sp.pres,value = .points.crs)}
    if (sf::st_crs(data.sp.pres) != terra::crs(ras.hii)){
      if (is.na (terra::crs (ras.hii)) ) {sf::st_set_crs (data.sp.pres,value = NA)  ; warning ('ASSUMING Points and HumanRaster data have the same projection')}
      if (!is.na (terra::crs (ras.hii)) ) {data.sp.pres <- sf::st_transform(data.sp.pres, terra::crs(ras.hii) )}
    }

    hii.sp.pres <-terra::extract (ras.hii, y = data.sp.pres, cells=FALSE, xy=TRUE)
    row.id.hii.NA <- which(is.na(hii.sp.pres[,2]))
    if (length (row.id.hii.NA) != 0) {
      output.comments <- rep (NA,length =nrow (xydat))
      output.comments [row.id.hii.NA] <- paste0('No human filter available for this record;' )
    }
    
    hii.sp.value = hii.sp.pres[,2]

    hii.sp.pres <-( hii.sp.value >= .th.human.influence)

    out$humanDetection_HumanInfluence_value=hii.sp.value
    out$humanDetection_HumanInfluence_test=hii.sp.pres
    out$humanDetection_HumanInfluence_comments= paste('Threshold influence=',.th.human.influence)

  }

  #start urban areas
  if (any (method %in% c('urban','all')))  {

    #check ref exists
    newoutdir = paste0(tempdir(),'/spatialData')
    alreadyDownloaded = file.exists (x = paste0(newoutdir,'/NE_urbanareas.shp'))
    if (alreadyDownloaded) myRef = sf_object <- st_read(dsn = newoutdir,layer = 'NE_urbanareas')
    if (!alreadyDownloaded) myRef= NULL 
    cc_urb_test =  .cc_urb_occTest(x = df, lon =xf,lat=yf ,value='flagged', verbose = FALSE,ref = myRef,outdir = output.dir )
    cc_urb_test <- (!  cc_urb_test) * 1
    out$humanDetection_UrbanAreas_value <- cc_urb_test
    out$humanDetection_UrbanAreas_test <- as.logical(cc_urb_test) 
    out$humanDetection_UrbanAreas_comments <- c('Urban areas from rnaturalearth')
  }
  
  #compute score
  out$humanDetection_score <-  .gimme.score (out)


  return (out)
}

# Biodiversity institutions ====
#' @title Detect biodiversity institutions  
#' @description Detect occurrences potentially in biodiversity institutions using different methods
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param lf The locality field in df.
#' @param method charcter. Vector of methods to be used. See details. Default 'all'
#' @param .points.crs crs argument for df
#' @param do logical. Should range analysis be performed?
#' @param verbose logical. Print messages? Default FALSE
#' @return data.frame
#' @details current implemented methods are : \cr
#' "fromCoordinates" (use information on localities of institutions, from \link[CoordinateCleaner]{CoordinateCleaner-package}) \cr
#' "fromBotanicLocalityName" (using information on locality names that include keywords of institutions)
#' @keywords internal 
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @seealso \link[CoordinateCleaner]{cc_inst} \link[CoordinateCleaner]{cc_gbif} 
#' @export
institutionLocality <- function (df,
                                 xf,
                                 yf,
                                 lf,
                                 method='all',
                                 .points.crs,
                                 do=TRUE, verbose=FALSE
                                 ){

  #output table
  out = data.frame (institutionLocality_fromBotanicLocalityName_value=NA,
                    institutionLocality_fromBotanicLocalityName_test=NA,
                    institutionLocality_fromBotanicLocalityName_comments=NA,
                    institutionLocality_fromCoordinates_value=NA,
                    institutionLocality_fromCoordinates_test=NA,
                    institutionLocality_fromCoordinates_comments=NA,
                    institutionLocality_score=NA
                    )[1:nrow (df),]

  row.names (out) <- NULL

  if (!do) {return (out)}
  
  #method fromBotanicLocalityName
  if (any (method %in% c('fromBotanicLocalityName','all'))) {

    #if no locality data is provided
    if(is.null (lf)) {
      if(verbose) print ('No locality in input data')
    }

    if(!is.null (lf)) {
      #if locality data is  provided
      loc=df[,lf]
      potential.names <- c('botanic', 'botanische', 'botanico', 'jardin', 'garden', 'botanical')
      localities <- as.character (loc)
      bot.garden <-sapply (localities, function (x) {
        if (is.na(x)){return (NA)}
        if (x==''){return (0)}
        splits = c(' ',',',';',':','/','\\\\','[.]')
        if (x %in% splits ) {return (0)}

        m <- tolower(.multiple.strsplit(x,multiple.splits = splits))
        potential.bot.garden <- sum ((m %in% potential.names) *1)
        potential.bot.garden <- ifelse(potential.bot.garden>0, 1, 0)
        return (potential.bot.garden)
        
      })
      bot.garden <- as.numeric (bot.garden)
      
      out$institutionLocality_fromBotanicLocalityName_value = bot.garden
      out$institutionLocality_fromBotanicLocalityName_test = as.logical (bot.garden)
      out$institutionLocality_fromBotanicLocalityName_comments = 'likely botanical garden'
      
    }
    

  }
  
  #method from Coordinates
  if (any (method %in% c('fromCoordinates','all'))) {
    cc_inst_test = CoordinateCleaner::cc_inst(x = df, lon =xf,lat=yf ,value='flagged', verbose = FALSE )
    cc_inst_test <- (!  cc_inst_test) * 1
    cc_gbif_test = CoordinateCleaner::cc_gbif(x = df, lon =xf,lat=yf ,value='flagged', verbose = FALSE )
    cc_gbif_test <- (!  cc_gbif_test) * 1

    cc_gbif_comments = ifelse(cc_gbif_test==1,'GBIFheadquarters',NA)
    cc_inst_comments = ifelse(cc_inst_test==1,'BiodivInst',NA)
    
    out$institutionLocality_fromCoordinates_value = ifelse ((cc_gbif_test==1 | cc_inst_test==1),1,0)
    out$institutionLocality_fromCoordinates_test = as.logical (ifelse ((cc_gbif_test==1 | cc_inst_test==1),1,0))
    out$institutionLocality_fromCoordinates_comments =  .paste3(cc_gbif_comments,cc_inst_comments)

  }

  out$institutionLocality_score <-  .gimme.score (out)

  return (out)

}


# geoOutliers ====
#' @title Detect geographic outliers
#' @description Detect geographical outliers using several tests
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param .alpha.parameter parameter setting for alphahull
#' @param .distance.parameter numeric. Maximum distance allowed. Default to 1000
#' @param .medianDeviation.parameter  numeric. Deviation parameter to . Default to 0.1
#' @param method charcter. Vector of methods to be used. See details. Default 'all'
#' @param .points.crs crs character. Indicate coordinate reference system
#' @param do logical. Should tests be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @details Methods implented are 'alphaHull', detecting species outside an alphaHull level (default 2), \cr
#' 'alphaHull' \cr
#' 'distance' corresponds to 'distance' method in  CoordinateCleaner::cc_outl(method='distance'). See ?CoordinateCleaner::cc_outl \cr
#' 'median' corresponds to 'mad' method in  CoordinateCleaner::cc_outl(method='mad'). See ?CoordinateCleaner::cc_outl \cr
#' 'quantSamplingCorrected' corresponds to 'mad' method in  CoordinateCleaner::cc_outl(method='quantile'). See ?CoordinateCleaner::cc_outl \cr
#' 'grubbs' implements Grubbs test to find spatial outliers. See ?findSpatialOutliers for details \cr
#' @return data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr), Cory Merow (cmerow@@gmail.com)
#' @seealso getPointsOutAlphaHull, \link[CoordinateCleaner]{cc_outl}, findSpatialOutliers
#' @export

geoOutliers         <- function (df,
                                xf,
                                yf,
                                .alpha.parameter,
                                .distance.parameter=1000,
                                .medianDeviation.parameter=5,
                                .samplingIntensThreshold.parameter=0.1,
                                .mcp_percSample =95,
                                method = 'all',
                                .points.crs ,
                                do=TRUE, verbose=FALSE){
  
 
  #build dataframe of all potential results
  out = data.frame (geoOutliers_alphaHull_value=NA,
                    geoOutliers_alphaHull_test=NA,
                    geoOutliers_alphaHull_comments = NA,
                    
                    geoOutliers_distance_value=NA,
                    geoOutliers_distance_test=NA,
                    geoOutliers_distance_comments = NA,
                    
                    geoOutliers_median_value=NA,
                    geoOutliers_median_test=NA,
                    geoOutliers_median_comments=NA,
                    
                    geoOutliers_quantileSamplingCorr_value=NA,
                    geoOutliers_quantileSamplingCorr_test=NA,
                    geoOutliers_quantileSamplingCorr_comments=NA,
                    
                    geoOutliers_Grubbs_value=NA,
                    geoOutliers_Grubbs_test=NA,
                    geoOutliers_Grubbs_comments=NA,
                    
                    geoOutliers_Dixon_value=NA,
                    geoOutliers_Dixon_test=NA,
                    geoOutliers_Dixon_comments=NA,
                    
                    geoOutliers_Rosner_value=NA,
                    geoOutliers_Rosner_test=NA,
                    geoOutliers_Rosner_comments=NA,

                    geoOutliers_mcp_value=NA,
                    geoOutliers_mcp_test=NA,
                    geoOutliers_mcp_comments=NA,
                    
                    geoOutliers_score=NA
                    
                    
                    )[1:nrow (df),]

  row.names (out) <- NULL

  if (!do) {return (out)}

  #prepare df for different methods of outliers in different packages
  df$species <- 'MyFakeSp'
  rownames(df) <-1:nrow(df)
  xydat <- df[,c(xf,yf)]
  if (nrow (xydat)<5) {warning('geoOutliers not performed N<5'); return(out)}

  if (any (method %in% c('alphaHull','all'))){

    #We typically will consider by default an alpha 2 -like ALA
    if (nrow (df) > 5  & !is.na(.alpha.parameter)) {
      points.outside.alphahull <- try (getPointsOutAlphaHull  (xydat,alpha = .alpha.parameter), silent=TRUE)
      if (inherits(points.outside.alphahull, 'try-error')) points.outside.alphahull <- rep (NA,length.out=nrow (xydat))
      out.comments <- paste0('GeoIndicator alphaHull (Alpha=',.alpha.parameter,')')
    }


    if (nrow (df) <=5 ) {
      points.outside.alphahull <- rep (NA, nrow(df) )
      out.comments <- paste0('N<=5. Analysis not run')
    }
    
    #write results into the results dataframe
    out$geoOutliers_alphaHull_value <- points.outside.alphahull
    out$geoOutliers_alphaHull_test <- as.logical (points.outside.alphahull)
    out$geoOutliers_alphaHull_comments <- out.comments

    ### NOT IMPLEMENTED YET
    #This is a bit experimental, when the user does not provide a predetermined value of alpha
    # we perform some analysis to get an "optimal" alpha parameter.
    # not implemented yet
    # if (nrow (df) >5 & is.na(.alpha.parameter) ) {
    #
    #   #build alpha hulls and determine alpha parameter for each from 1 to 10
    #   ah.hulls.list <- lapply (1:10, function (a){
    #
    #     try(ah <- alphahull::ahull (x = df[,xf], y= df[,yf], alpha=a), silent = TRUE)
    #     if (exists("ah")){return (ah)}
    #
    #   })
    #
    #   ah.hulls.thresholds <- lapply (1:10, function (i) {
    #     if((is.null(ah.hulls.list[[i]]) == FALSE) & (is.list (ah.hulls.list[[i]]) == TRUE)  ) {inahull(ah.hulls.list[[i]] , as.matrix(xydat) ) * i}
    #   })
    #
    #   is.not.null <- lapply (ah.hulls.thresholds, function (x) !is.null(x))
    #   ah.hulls.thresholds <- ah.hulls.thresholds[unlist(is.not.null)]
    #   ah.hulls.thresholds <- Reduce (cbind, ah.hulls.thresholds)
    #
    #   #check if error in alphahulls, if not continue
    #   if ( is.null (ah.hulls.thresholds)) {
    #     ah.hulls.level <- rep (NA, nrow(df) )
    #     for (i in 1:nrow(df)) {df$comments [i] <- paste0(df$comments [i],'Error calculating AlphaHull;' )}
    #
    #   }
    #   if (!is.null (ah.hulls.thresholds)) {
    #     ah.hulls.level <- apply (ah.hulls.thresholds, 1, function (x){
    #
    #       x [which (x==0)] <- NA
    #       m <- min (x, na.rm=TRUE)
    #       out <- ifelse (m==Inf, NA, m)
    #       return (out)
    #
    #
    #     })
    #     ah.hulls.level[which(is.na(ah.hulls.level))] <- 9999
    #
    #     #build dataframe for optimal alpha value choice
    #     ah.dataframe <- lapply (1:length(ah.hulls.list), function (i){
    #
    #       ah <- ah.hulls.list[[i]]
    #       if (exists("ah") & is.null(ah) == FALSE & is.list(ah)==TRUE ){
    #         perc.points.in.hull <- round ((sum (alphahull::inahull(ah,as.matrix(xydat)))/nrow (xydat))*100,digits=2)
    #         area.hull <-alphahull::areaahull(ah)
    #         return (data.frame (perc.points.in.hull, area.hull))
    #
    #       } else {
    #         perc.points.in.hull <- NA
    #         area.hull <-   NA
    #         return (data.frame (perc.points.in.hull, area.hull))
    #
    #       }
    #
    #
    #
    #
    #     })
    #     ah.dataframe <- Reduce (rbind,ah.dataframe)
    #     ah.dataframe$alpha <- 1:10
    #     ah.dataframe$optim.alpha.values <- ah.dataframe$perc.points.in.hull/ah.dataframe$area.hull
    #     dat.ts <- stats::ts(ah.dataframe$area.hull, frequency=1)
    #     ah.dataframe$A.increase <- c (dat.ts/stats::lag(x = dat.ts, 1) - 1,NA)
    #
    #     #choose what I would call the optimal alpha based on two criteria
    #
    #     alpha.choice <- subset (ah.dataframe, subset=perc.points.in.hull > 90)
    #     if (nrow (alpha.choice) == 0) {print('No optimal alpha choice with the criteria used. Defaulting to alpha = 2')
    #       .alpha.parameter <- 2
    #
    #     }
    #
    #     if (nrow (alpha.choice) >= 0){
    #       alpha.choice <- alpha.choice[order(alpha.choice$optim.alpha.values,decreasing = TRUE,na.last = TRUE),]
    #       alpha.choice <- alpha.choice[1,'alpha']
    #       if(verbose) print (paste("Optimal alpha selected =", alpha.choice, "for alpha geoindicator of outliers"))
    #
    #       .alpha.parameter <- alpha.choice
    #
    #     }
    #
    #
    #   }
    #
    #
    #
    #
    # }
    #
    #



  }

  if (any (method %in% c('distance','all'))){
    #create fake species list for Coordinate cleaner
    df$species <- 'MyFakeSp'
    cc_dist_test =  .cc_outl_occTest(x = df,lon = xf,lat = yf,method = 'distance',value = 'flagged',tdi = .distance.parameter, verbose=FALSE)
    cc_dist_test <- (!  cc_dist_test) 
    out$geoOutliers_distance_value  = cc_dist_test * 1
    out$geoOutliers_distance_test  = cc_dist_test
    out$geoOutliers_distance_comments = rep(paste('dist.threshold=',.distance.parameter,'m'),times=nrow(out))
    if (nrow (df)<7) {
      out$geoOutliers_distance_value = rep(NA,times=nrow(out))
      out$geoOutliers_distance_test = rep(NA,times=nrow(out))
      out$geoOutliers_distance_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))
    }
  }
  
  if (any (method %in% c('median','all'))){

    cc_med_test =  .cc_outl_occTest(x = df,lon = xf,lat = yf,method = 'mad',value = 'flagged',mltpl =  .medianDeviation.parameter , verbose=FALSE)
    cc_med_test <- (!  cc_med_test) * 1
    out$geoOutliers_median_value = cc_med_test
    out$geoOutliers_median_test  = as.logical(cc_med_test)
    out$geoOutliers_median_comments = rep(paste('medDeviation=', .medianDeviation.parameter),times=nrow(out))

    if (nrow (df)<7) {
      out$geoOutliers_median_value = rep(NA,times=nrow(out))
      out$geoOutliers_median_test = rep(NA,times=nrow(out))
      out$geoOutliers_median_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))
    }
  }

  if (any (method %in% c('quantSamplingCorrected','quantileSamplingCorr','all'))){
    cc_qsc_tst =  .cc_outl_occTest(x = df,lon = xf,lat = yf,method='quantile',value = 'flagged',
                                            mltpl =  .medianDeviation.parameter,
                                            sampling_thresh = .samplingIntensThreshold.parameter,
                                             verbose=FALSE)
    cc_qsc_tst <- (!  cc_qsc_tst) * 1
    out$geoOutliers_quantileSamplingCorr_value = cc_qsc_tst
    out$geoOutliers_quantileSamplingCorr_test = as.logical(cc_qsc_tst)
    out$geoOutliers_quantileSamplingCorr_comments = rep(paste('IQRmultiplier=', .medianDeviation.parameter,';SamplIntensityThres=',.samplingIntensThreshold.parameter),times=nrow(out))
    if (nrow (df)<7) {
      out$geoOutliers_quantileSamplingCorr_value = rep(NA,times=nrow(out))
      out$geoOutliers_quantileSamplingCorr_test = rep(NA,times=nrow(out))
      out$geoOutliers_quantileSamplingCorr_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))
      }
  }

  if (any (method %in% c('grubbs','all','Grubbs'))) {
    spdf = sf::st_as_sf(xydat,coords =c(xf,yf),crs =.points.crs )
    if (nrow (spdf)<=5) {
      out$geoOutliers_Grubbs_value = rep(NA,times=nrow(out))
      out$geoOutliers_Grubbs_test = rep(NA,times=nrow(out))
      out$geoOutliers_Grubbs_comments = rep(paste('N<5. Analysis not run'),times=nrow(out))
      }
    if (nrow (spdf)>5){
      grubbs_tst = rep(0,times=nrow(out))
      grubbs.outl = try(findOutlyingPoints(pres = spdf,spOutliers = T,envOutliers = F,verbose = F),silent=TRUE)
      if (inherits(grubbs.outl,what = 'try-error')) {
        grubbs_tst = rep(NA,times=nrow(out)) 
        grubbs_comment <- rep (paste('Error in findSpatialOutliers'),times=nrow(out))
      }
      if (!inherits(grubbs.outl,what = 'try-error')) {
        grubbs_tst <- as.numeric (grubbs.outl$spOutlier)
        grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))
      }
      out$geoOutliers_Grubbs_value = grubbs_tst
      out$geoOutliers_Grubbs_test = as.logical(grubbs_tst)
      out$geoOutliers_Grubbs_comments = grubbs_comment
    }
  }
  
  if (any (method %in% c('dixon','all','Dixon'))) {
    spdf = sf::st_as_sf(xydat,coords =c(xf,yf),crs =.points.crs )
    if (nrow (spdf)<=5) {
      out$geoOutliers_Dixon_value = rep(NA,times=nrow(out))
      out$geoOutliers_Dixon_test = rep(NA,times=nrow(out))
      out$geoOutliers_Dixon_comments = rep(paste('N<5. Analysis not run'),times=nrow(out))
    }
    if (nrow (spdf)>5){
      dixon_tst = rep(0,times=nrow(out))
      dixon.outl = try(findOutlyingPoints(pres = spdf,spOutliers = T,method = 'dixon',envOutliers = F,verbose = F,),silent=TRUE)
      if (inherits(dixon.outl,what = 'try-error')) {
        dixon_tst = rep(NA,times=nrow(out)) 
        dixon_comment <- rep (paste('Error in findSpatialOutliers'),times=nrow(out))
      }
      if (!inherits(dixon.outl,what = 'try-error')) {dixon_tst <- as.numeric (dixon.outl$spOutlier)}
      out$geoOutliers_Dixon_value = dixon_tst
      out$geoOutliers_Dixon_test = as.logical(dixon_tst)
      out$geoOutliers_Dixon_comments = NA
    }
  }
  
  if (any (method %in% c('rosner','all','Rosner'))) {
    spdf = sf::st_as_sf(xydat,coords =c(xf,yf),crs =.points.crs )
    if (nrow (spdf)<=5) {
      out$geoOutliers_Rosner_value = rep(NA,times=nrow(out))
      out$geoOutliers_Rosner_test = rep(NA,times=nrow(out))
      out$geoOutliers_Rosner_comments = rep(paste('N<5. Analysis not run'),times=nrow(out))
    }
    if (nrow (spdf)>5){
      rosner_tst = rep(0,times=nrow(out))
      rosner.outl = try(findOutlyingPoints(pres = spdf,spOutliers = T,method = 'rosner',
                                           envOutliers = F,verbose = F,),silent=TRUE)
      if (inherits(rosner.outl,what = 'try-error')) {
        rosner_tst = rep(NA,times=nrow(out)) 
        rosner_comment <- rep (paste('Error in findSpatialOutliers'),times=nrow(out))
      }
      if (!inherits(rosner.outl,what = 'try-error')) {rosner_tst <- as.numeric (rosner.outl$spOutlier)}
      out$geoOutliers_Rosner_value = rosner_tst
      out$geoOutliers_Rosner_test = as.logical(rosner_tst)
      out$geoOutliers_Rosner_comments = NA
    }
  }
  
  if (any (method %in% c('mcp','all'))) {
      xydat_mcpTest <- mcp_test(xydat,xf,yf,percentage=.mcp_percSample)
      out$geoOutliers_mcp_value = xydat_mcpTest
      out$geoOutliers_mcp_test = as.logical(xydat_mcpTest)
      if (all(is.na(out$geoOutliers_mcp_value))) {
        out$geoOutliers_mcp_comments = rep(paste('N<5. Analysis not run'),times=nrow(out)) 
      }
    }
    
  out$geoOutliers_score <-  .gimme.score (out)

  return (out)
}

# env outliers ====
#' @title Detect environmental outliers
#' @description Detect environmental outliers using different methods
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param .r.env raster. Raster with environmental data
#' @param .th.perc.outenv numeric. Value from 0 to 1 to identify the rate of variables not passing the test to consider the record an environmental outlier.
#' @param .sp.name Species name
#' @param method charcter. Vector of methods to be used. See details. Default 'all'
#' @param .points.crs crs character. Indicate coordinate reference system
#' @param do logical. Should tests be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @return data.frame
#' @details Implemented methods are :\cr
#' 'bxp' based on boxplot distribution (1.5 Interquartile range) in a variable to identify the record as outlier (\link[biogeo]{outliers}) \cr
#' 'Grubbs' based on Grubbs test. See ?findEnvOutliers \cr
#' Regardless of the method, the function already tests for missing environmental variables and it outputs the result in the output data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr), Cory Merow (cmerow@@gmail.com). Mark Robertson (biogeo package)
#' @seealso findEnvOutliers \link[biogeo]{outliers}
#' @import stats 
#' @export
envOutliers  <- function (
                          df,
                          xf,
                          yf,
                          .r.env,
                          .th.perc.outenv ,
                          .sp.name ,
                          method='all',
                          .points.crs ,
                          do=TRUE, verbose=FALSE){

  #output results
  out = data.frame (envOutliers_missingEnv_value = NA,
                    envOutliers_missingEnv_test = NA,
                    
                    envOutliers_bxp_value = NA,
                    envOutliers_bxp_test = NA,
                    envOutliers_bxp_comments = NA,
                    
                    envOutliers_Grubbs_value=NA,
                    envOutliers_Grubbs_test=NA,
                    envOutliers_Grubbs_comments=NA,
                    
                    envOutliers_Dixon_value=NA,
                    envOutliers_Dixon_test=NA,
                    envOutliers_Dixon_comments=NA,
                    
                    envOutliers_Rosner_value=NA,
                    envOutliers_Rosner_test=NA,
                    envOutliers_Rosner_comments=NA,
                    
                    
                    envOutliers_jack_value=NA,
                    envOutliers_jack_test=NA,
                    envOutliers_jack_comments=NA,
                    
                    envOutliers_svm_value=NA,
                    envOutliers_svm_test=NA,
                    envOutliers_svm_comments=NA,
                    
                    envOutliers_rf_value=NA,
                    envOutliers_rf_test=NA,
                    envOutliers_rf_comments=NA,
                    
                    envOutliers_rfout_value=NA,
                    envOutliers_rfout_test=NA,
                    envOutliers_rfout_comments=NA,
                    
                    envOutliers_lof_value=NA,
                    envOutliers_lof_test=NA,
                    envOutliers_lof_comments=NA,
                    
                    envOutliers_score=NA
                    )[1:nrow (df),]

  if (!do) {return (out)}
  #prepare df for different methods of outliers in different packages
  df$species <- 'MyFakeSp'
  rownames(df) <-1:nrow(df)
  xydat <- df[,c(xf,yf)]
  if (nrow (xydat)<5) {
    warning('envOutliers not performed N<5')
    
    return(out)}
  dat.environment <- terra::extract(x = .r.env, y= df[,c(xf,yf)], xy = T,ID=F)
  #subs coord of the cell by coord of the point
  dat.environment$x <- df[,xf]
  dat.environment$y <- df[,yf]
  missingEnvironment <- (!stats::complete.cases(dat.environment))
  out$envOutliers_missingEnv_value <- missingEnvironment * 1
  out$envOutliers_missingEnv_test <- missingEnvironment
  out$envOutliers_missingEnv_comment <- ''
  
  if (any (missingEnvironment)){
    ids = which (missingEnvironment==TRUE)
    for (i in ids){
      de = dat.environment[i,]
      layersWithNA = which (is.na(de))
      out$envOutliers_missingEnv_comment[i] <- paste ("layers with NA:",paste(names (de)[layersWithNA],collapse =','))
    }
  }
  if (any (method %in% c('bxp','all'))) {

    outlier.method.biogeo = '-bxp'
    # dat.environment <- extract(x = .r.env, y= df[,c(xf,yf)], df = TRUE)
    # missingEnvironment <- (!stats::complete.cases(dat.environment))*1
    species <- rep(.sp.name, nrow(dat.environment))
    dups <- rep(0, nrow(dat.environment))
    environmental.columns <- names (.r.env) #2:ncol(dat.environment)
    a.df <- lapply (environmental.columns, function (i){

      ev<-dat.environment[,i]

      a<- try (biogeo::outliers(rid=1:nrow(dat.environment), species, dups, ev), silent=TRUE)
      if (exists('a') & all(class(a)!='try-error')){
        a<- as.data.frame (a)
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}

      if (!exists('a') | any(class(a)=='try-error')){
        #in case there is no variation across the variable and then outliers function throws an error
        a <- data.frame (col1=rep (0,nrow(dat.environment)), col2 = rep (0,nrow(dat.environment)))
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}
    })
    a.df <- do.call (cbind, a.df)
    choice.of.method.df <- a.df [,grep (names(a.df),pattern = outlier.method.biogeo)]
    if (inherits(.r.env,"SpatRaster")) {
      if (terra::nlyr(.r.env) == 1) outlier.level <- choice.of.method.df
      if (terra::nlyr(.r.env) > 1)  outlier.level <- round (rowMeans(choice.of.method.df) ,digits=2)
      }
    outlier.env <- (outlier.level >= .th.perc.outenv)*1
    comments.outlier.env <- paste (paste("Outlier",outlier.method.biogeo,".level="), outlier.level)
    out$envOutliers_bxp_value = outlier.env
    out$envOutliers_bxp_test = as.logical(outlier.env)
    out$envOutliers_bxp_comments = comments.outlier.env
  }
  if (any (method %in% c('grubbs','all','Grubbs')))  {
    spdf <- sf::st_as_sf(dat.environment,coords=c('x','y'),crs = .points.crs)
    env.grubbs_tst = rep(0,times=nrow(out))
    env.grubbs.outl = try(findOutlyingPoints(pres = spdf,
                                             spOutliers = F,envOutliers = T,
                                             method = 'grubbs',verbose = F),silent=TRUE)
    
    if (inherits(env.grubbs.outl,what = 'try-error')) {
      env.grubbs_tst = rep(NA,times=nrow(out)) 
      env.grubbs_comment <- rep (paste('Error in findSpatialOutliers'),times=nrow(out))
    }
    
    if (!inherits(env.grubbs.outl,what = 'try-error')) {
      env.grubbs_tst <- as.numeric (env.grubbs.outl$envOutlier)
      #grubbs_tst [grubbs.outl] = 1
      env.grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))
    }
    
    out$envOutliers_Grubbs_value = env.grubbs_tst * ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1)
    out$envOutliers_Grubbs_test = as.logical(env.grubbs_tst * ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1))
    out$envOutliers_Grubbs_comments = env.grubbs_comment

  }
  if (any (method %in% c('dixon','all','Dixon')))  {
    spdf <- sf::st_as_sf(dat.environment,coords=c('x','y'),crs = .points.crs)
    env.dixon.outl = try(findOutlyingPoints(pres = spdf,
                                             spOutliers = F,envOutliers = T,
                                             method = 'dixon',verbose = F),silent=TRUE)
    
    if (inherits(env.dixon.outl,what = 'try-error')) {
      env.dixon_tst = rep(NA,times=nrow(out)) 
      env.dixon_comment <- rep (paste('Error in findOutlyingPoints'),times=nrow(out))
    }
    if (!inherits(env.dixon.outl,what = 'try-error')) {
      env.dixon_tst <- as.numeric (env.dixon.outl$envOutlier)
      env.dixon_comment <- rep (NA,times=nrow(out))
    }
    out$envOutliers_Dixon_value = env.dixon_tst * ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1)
    out$envOutliers_Dixon_test = as.logical(env.dixon_tst * ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1))
    out$envOutliers_Dixon_comments = env.dixon_comment
    
  }
  if (any (method %in% c('rosner','all','Rosner')))  {
    spdf <- sf::st_as_sf(dat.environment,coords=c('x','y'),crs = .points.crs)
    env.rosner_tst = rep(0,times=nrow(out))
    env.rosner.outl = try(findOutlyingPoints(pres = spdf,
                                            spOutliers = F,envOutliers = T,
                                            method = 'rosner',verbose = F),silent=TRUE)
    if (inherits(env.rosner.outl,what = 'try-error')) {
      env.rosner_tst = rep(NA,times=nrow(out)) 
      env.rosner_comment <- rep ('Error in findOutliers Rosner',times=nrow(out))
    }
    if (!inherits(env.rosner.outl,what = 'try-error')) {
      env.rosner_tst <- as.numeric (env.rosner.outl$envOutlier)
      env.rosner_comment <- rep (NA,times=nrow(out))
    }
    out$envOutliers_Rosner_value = env.rosner_tst * ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1)
    out$envOutliers_Rosner_test = as.logical(env.rosner_tst* ifelse (out$envOutliers_missingEnv_test,yes = NA,no = 1))
    out$envOutliers_Rosner_comments = env.rosner_comment
    
  }
  #compute if any flexsdm methods
  if (any (method %in% c('jack','all','svm','rf','rfout','lof'))){
    dat.environment_flexsdm_pres <- dat.environment %>% 
      tibble::add_column(pr_ab=1)
    #simulate absences 
    xy_vect <- terra::vect(dat.environment, geom=c("x", "y"),crs=terra::crs("epsg:4326"))
    d <- median (terra::distance(xy_vect),na.rm=T)
    comment_dist_flexsdm <- paste('buffer distance',d,'m')
    b_xy <- terra::buffer(xy_vect,width=d)
    r.env_crop <- terra::crop(.r.env,terra::ext(b_xy))
    rr_bxy <- terra::mask (x = r.env_crop,mask = b_xy)
    cid_pres <- terra::cells(rr_bxy,xy_vect)
    pres_mask  =terra::rasterize(as.matrix(dat.environment_flexsdm_pres[,c('x','y')]),
                                 rr_bxy,
                                 vals=1,
                                 background=0)
    pres_mask = terra::subst(pres_mask,
                             from=c(1,0),
                             to=c(NA,1))
    rr_bxy_presmask <- pres_mask * rr_bxy
    names (rr_bxy_presmask) <- names (.r.env)
    rrCellsNotNA = min (terra::global(rr_bxy_presmask,'notNA'))
    Nabs_size = min (length(cid_pres)*3,
                 2000,
                 rrCellsNotNA-length(cid_pres),
                 rrCellsNotNA*0.6)
    #check if you have the minimum N of abs to compute abs methods
    if (Nabs_size *1.2 >= length(cid_pres)){
      dat.environment_flexsdm_abs = terra::spatSample (rr_bxy_presmask,
                                                       values=T,
                                                       as.df=T,
                                                       xy=T,
                                                       na.rm=T,
                                                       replace=T,
                                                       method='random',
                                                       size=Nabs_size)
      dat.environment_flexsdm_abs$pr_ab <- 0
      dat.environment_flexsdm <- dplyr::bind_rows (dat.environment_flexsdm_pres,
                                                   dat.environment_flexsdm_abs) 
      dat.environment_flexsdm$id <- 1:nrow (dat.environment_flexsdm)
      comment_abs_flexsdm <- paste('N absences =',Nabs_size)
    }
    if (Nabs_size *1.2 < length(cid_pres)){
      dat.environment_flexsdm <- dat.environment_flexsdm_pres 
      dat.environment_flexsdm$id <- 1:nrow (dat.environment_flexsdm)
      comment_abs_flexsdm <- 'PresAbs methods not computed. Too many presences'
    }
    #using local import and modf of flexsdm
    flexsdm_out <- env_outliers(dat.environment_flexsdm,x = 'x',y = 'y',pr_ab = 'pr_ab',id='id',env_layer = .r.env)
  }
  #add flexsdm outputs
  if (any (method %in% c('jack','all')))  {
    envOutliers_jack_value <- flexsdm_out %>% 
      dplyr::filter (pr_ab==1) %>% 
      dplyr::pull(dplyr::ends_with('_jack'))
    envOutliers_jack_test <- as.logical (envOutliers_jack_value)
  }
  if (any (method %in% c('svm','all')))  {
    envOutliers_svm_value <- flexsdm_out %>% 
      dplyr::filter (pr_ab==1) %>%
      dplyr::pull(dplyr::ends_with('_svm'))
    envOutliers_svm_test <- as.logical (envOutliers_svm_value)
    envOutliers_svm_comment <- comment_abs_flexsdm
  }
  if (any (method %in% c('rf','all')))  {
    envOutliers_rf_value <- flexsdm_out %>%
      dplyr::filter (pr_ab==1) %>%
      dplyr::pull(dplyr::ends_with('_rf'))
    envOutliers_rf_test <- as.logical (envOutliers_rf_value)
    envOutliers_rf_comment <- comment_abs_flexsdm
  }
  if (any (method %in% c('rfout','all')))  {
    envOutliers_rfout_value <- flexsdm_out %>% 
      dplyr::filter (pr_ab==1) %>%
      dplyr::pull(dplyr::ends_with('_rfout'))
    envOutliers_rfout_test <- as.logical (envOutliers_rfout_value)
  }
  if (any (method %in% c('lof','all')))  {
    envOutliers_lof_value <- flexsdm_out %>% 
      dplyr::filter (pr_ab==1) %>%
      dplyr::pull(dplyr::ends_with('_lof'))
    envOutliers_lof_test <- as.logical (envOutliers_lof_value)
  }
  out$envOutliers_score <-  .gimme.score (out)
  return (out)
}

# geoEnvAccuracy ====
#' @title Coordinate accuracy
#' @description Detect records with low accuracy affecting environmental values
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param af character. column name in df containing the coordinate uncertainty value (in the same)
#' @param dsf character. column name in df containing the dataset to which the record belongs to (e.g. Forest Inventory of Spain)
#' @param ef character. column name in df containing the registered elevation for the record. 
#' @param method character. Vector of methods to be used. See details. Default 'all'
#' @param r.env raster. Raster with environmental data
#' @param accept.threshold.cell numeric. Acceptance threshold for how much percentage of the Area of uncertainty in the cell we want to accept. Default to 0.5
#' @param accept.threshold.env numeric. Default 0.5
#' @param bearing.classes numeric. Default to 10.
#' @param distance.classes integer. Default to 5.
#' @param env.quantiles numeric. Default to c(0.3,0.7)
#' @param elev.threshold numeric. Default to 100
#' @param raster.elevation numeric. Default to 100
#' @param do logical. Should tests be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @param doParallel logical. Should computation use parallel functions? Default FALSE
#' @param mc.cores numeric. How many cores to use? (used when doParallel = TRUE). Default 2 
#' @return data.frame
#' @keywords Analysis 
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr), A Zizka (CoordinateCleaner package)
#' @details Geoenvironmental accuracy function will implement differnt methods to assess occurrence accuracy in environmnental and geographic space.\cr
#' Current implmented methods are:
#' 'lattice' : tests for lattice arrangement in occurrence datasets. Borrowed from \link[CoordinateCleaner]{cd_round} . \cr
#' 'elevDiff' : assess the elevation difference between a given raster (or automatically downloaded fro SRTM), and the elevation recorded. If differences >elev.threshold then the record is considered as a low accuracy threshold\cr
#' 'noDate' : assess whether there is a date or timestamp information in the record. \cr 
#' 'noDateFormatKnown' : assess whether the information in the timestamp agrees with different formatting of Dates. \cr 
#' 'outDateRange' : (not implemented) assess whether the record is within a user specified time frame. \cr 
#' 'percDiffCell' : assess whether the record may be falling in a different raster cell given an information of coordinate accuracy. \cr 
#' 'envDeviation' : assess whether the climate in a given record may be outside of the interval 30th-70th (default values) for a given variable due to coordinate uncertainty. \cr 
#' @seealso \link[CoordinateCleaner]{cd_round}
#' @family analysis
#' @examples 
#' #see examples in vignetteXtra-occTest
#' #library
#' #install_github ('vignetteXtra-occTest')
#' @export
#' @importFrom elevatr get_elev_point
#' @importFrom stats na.omit

geoEnvAccuracy  <- function (df,
                             xf,
                             yf,
                             af,dsf,ef,

                             
                             method='all',
                             
                             r.env,
                             accept.threshold.cell =0.5,
                             accept.threshold.env = 0.5,
                             
                             bearing.classes=10,
                             distance.classes=5,
                             env.quantiles=c(0.3,0.7),
                             
                             elev.threshold = 100,
                             raster.elevation=NULL,
                             verbose=FALSE,
                             do=TRUE,
                             doParallel=FALSE,
                             mc.cores= 2
){
  
  #output results
  out = data.frame (geoenvLowAccuracy_lattice_value = NA,
                    geoenvLowAccuracy_lattice_test = NA,
                    geoenvLowAccuracy_lattice_comments = NA,
                    
                    geoenvLowAccuracy_percDiffCell_value = NA,
                    geoenvLowAccuracy_percDiffCell_test = NA,
                    geoenvLowAccuracy_percDiffCell_comments = NA,
                    
                    geoenvLowAccuracy_envDiff_value =NA,
                    geoenvLowAccuracy_envDiff_test =NA,
                    geoenvLowAccuracy_envDiff_comments =NA,
                    
                    geoenvLowAccuracy_elevDiff_value =NA,
                    geoenvLowAccuracy_elevDiff_test =NA,
                    geoenvLowAccuracy_elevDiff_comments =NA,
            
                    geoenvLowAccuracy_score=NA)[1:nrow (df),]
  row.names(x = out) <- NULL
  if (!do) {return (out)}
  #check if need parallel
  os =  .get_os()
  if (doParallel == TRUE){
    mymclapply <- switch( Sys.info()[['sysname']],
                        Windows = {parallelsugar::mclapply},
                        Linux   = {parallel::mclapply},
                        Darwin  = {parallel::mclapply})
  }
  if (doParallel==FALSE) {mymclapply <- lapply}
  if (nrow(df)<100) {mymclapply <- lapply}
  
  #start method lattice
  if (any(method %in% c('lattice','all'))){
    
    if(!is.null(dsf)) {
      xydat = df[,c(xf,yf,dsf)]
    }
    if(is.null(dsf)) {
      print ('No dataset field provided. Assuming all records are from the same dataset')
      xydat = df[,c(xf,yf)]
      xydat$dataset = 'TemporaryDatasetName'
      dsf <- 'dataset'
    }
    
    #there are a lot of problems with this function
    cd_round_test = try({.cc_round_occTest(x = xydat,lon = xf,lat = yf,ds=dsf,value='flagged2',verbose=FALSE,graphs=FALSE)},silent=TRUE)
    #cd_round_test = CoordinateCleaner::cd_round(x = xydat,lon = xf,lat = yf,ds=dsf,value='flagged',verbose=FALSE,graphs=FALSE)
    out$geoenvLowAccuracy_lattice_value = (!cd_round_test) * 1
    out$geoenvLowAccuracy_lattice_test = (!cd_round_test)
    out$geoenvLowAccuracy_lattice_comments = ('Based on default values of CoordinateCleaner::cd_round')
  }
  
  #start methods requiring elevation
  #method elevation
  if (any(method %in% c('elevDiff','all'))) {
    if (!is.null(ef)){
      
      if (is.null (raster.elevation)) {
        if (verbose) print ('elevation raster not provided, downloading from SRTM')
        do.srtm.download = TRUE
      } else {do.srtm.download=FALSE}
      #if no elevation raster provided, download it
      if (do.srtm.download){
        xydat = df[,c(xf,yf)]
        assumedProj = "EPSG:4326"
        df_elev_aws <- elevatr::get_elev_point(xydat, 
                                                prj = assumedProj, 
                                                src = "aws")
      }
      #determining elevation test
      if (inherits(raster.elevation,"SpatRaster")){
        elev.xy <- terra::extract (raster.elevation, df[,c(xf,yf)]) 
        elev.xy <-elev.xy[,2]
        } 
      if (do.srtm.download & exists('df_elev_aws')) {
        elev.xy <-  df_elev_aws$elevation
        }
      elev.difference <- abs (df[,ef] - elev.xy)
      out$geoenvLowAccuracy_elevDiff_value <- elev.difference
      out$geoenvLowAccuracy_elevDiff_test <- as.logical (ifelse (elev.difference >= elev.threshold, 1, 0))
      out$geoenvLowAccuracy_elevDiff_comments <- paste ('Elev difference:', elev.difference)
      #skipping elevation test
      if (!inherits(raster.elevation,"SpatRaster") & !do.srtm.download) {
        if(verbose) warning ("raster.elevation is not a SpatRaster; skipping elevDiff method")
      }
    }
  }
  
  # #start time methods
  # #has time stamp? 
  # if (any(method %in% c('noDate','all'))) {
  #   if (!is.null(tf)){
  #     vecTime = df[,tf]
  #     dfEmptyVals = data.frame (na=is.na(vecTime),
  #                               #null=is.null(as.character(vecTime)==NULL),
  #                               isEmpty = (as.character(vecTime)==''),
  #                               isEmptySpace = (as.character(vecTime)== '  '))
  #     
  #     vecTimeLogi = apply(dfEmptyVals,MARGIN = 1,FUN = function (x){any (x==TRUE)})
  #     out$geoenvLowAccuracy_noDate_value = vecTimeLogi * 1
  #     out$geoenvLowAccuracy_noDate_test = vecTimeLogi 
  #     out$geoenvLowAccuracy_noDate_comments = 'checked for NA,NULL, blankSpace and tab' 
  #   }
  # }
  # #has format stamp? 
  # if (any(method %in% c('noDateFormatKnown','all'))) {
  #   if (!is.null(tf)){
  #     
  #     warning ('method noDateFormatKnown under development. Sorry')
  #     # vecTime = df[,tf]
  #     # 
  #     # dfVectime = data.frame(vecTime[1])
  #     # 
  #     # a = identify_dates(dfVectime)
  #     # dfVectime = data.table::as.data.table(dfVectime)
  #     # k = find_and_transform_dates(dfVectime)
  #     # 
  #     # 
  #     # out$geoenvLowAccuracy_noDate_value = vecTimeLogi * 1
  #     # out$geoenvLowAccuracy_noDate_test = vecTimeLogi 
  #     # out$geoenvLowAccuracy_noDate_comments = 'checked for NA,NULL, blankSpace and tab' 
  #   }
  #   
  # }
  # #is within time range ?
  # if (any(method %in% c('outDateRange','all'))) {
  #   if (!is.null(tf)){
  #     warning ('method outDateRange under development. Sorry')
  #     vecTime = df[,tf]
  #     
  #     
  #     
  #   }
  #   
  # }
  # 
  
  #Need coordinate uncertainty analysis ? If not gimme score and leave
  if (!any (method %in% c('percDiffCell','envDeviation','all')) ) {
    #write final score
    out$geoenvLowAccuracy_score <- .gimme.score (out)
    return (out)}
  #Need coordinate uncertainty analysis ? Can you do it?
  if (is.null(af) | all(is.na(df[,af]))) {
    print ('no coordinate accuracy/uncertainty field provided')
    #write final score
    out$geoenvLowAccuracy_score <- .gimme.score (out)
    return (out)
  }
  #Start methods requiring coordinate uncertainty
  if(length(af)>1) {df$new_accuracy = pmax(df[,af[1]],df[,af[2]],na.rm=TRUE); af = 'new_accuracy'}
  xydat = df[,c(xf,yf,af)]
  xydat$occCellids = terra::cellFromXY(r.env,as.data.frame (df[,c(xf,yf)]))
  #get cell IDs of buffer points
  cellIds.directions = mymclapply (1:nrow (xydat), function (i){
    distance.accuracy = as.numeric (xydat[i,c(af)])
    if (is.na (distance.accuracy) |is.na (xydat$occCellids[i])  ) {return (NA)}
    if (distance.accuracy>1){
      sequence.pts = seq(from=0, to=as.numeric(distance.accuracy) ,length.out = distance.classes +1)
      sequence.pts = sequence.pts [-1]
    } else {sequence.pts <- distance.accuracy}
    
    myAngles = seq(from=1,to=360,length.out=bearing.classes)
    out.bearing = sapply(myAngles, function (bearing){
      
      points.buffer.id = sapply(sequence.pts, function (dpoints){
        dest.point = geosphere::destPoint(p = xydat[i,c(xf,yf)],b=bearing,d=dpoints)
        dest.point.cellID = terra::cellFromXY(object = r.env, dest.point)
        #diff.2.target.cell = (xydat$occCellids[i] != dest.point.cellID)*1
        #diff.2.target.cell
        dest.point.cellID
      })
      #perc.points.buffer = mean (points.buffer,na.rm=TRUE)
      points.buffer.id
    })
    
    as.vector (out.bearing)
  })
  
  #start method percentage of Different cells around uncertainty
  if (any(method %in% c('percDiffCell','all'))){
    perc.diff <- sapply (1:nrow(xydat), function (i){
      pdiff <- mean ((cellIds.directions[[i]]==xydat$occCellids[i]) * 1, na.rm = TRUE)
      pdiff
    })
    out$geoenvLowAccuracy_percDiffCell_value = (perc.diff)
    out$geoenvLowAccuracy_percDiffCell_test = (perc.diff>accept.threshold.cell)
    out$geoenvLowAccuracy_percDiffCell_comments = paste('Nbufferlocations=',bearing.classes*distance.classes)
  }
  #start method percentage of environmental differences around uncertainty
  if (any(method %in% c('envDeviation','all'))) {
    
    uniqueCellBuff <- unique (unlist (cellIds.directions))
    uniqueCellBuff <- na.omit(uniqueCellBuff)
    if (terra::nlyr(r.env)>1) df.CellBuff = terra::extract(r.env,uniqueCellBuff)
    if (terra::nlyr(r.env)==1) {
      df.CellBuff = terra::extract(r.env,uniqueCellBuff)
      df.CellBuff = data.frame (ID=1:length(df.CellBuff),LAYER = df.CellBuff)
      names (df.CellBuff) = c('ID',names(r.env))
    } 
    df.CellBuff$cellID = uniqueCellBuff
    
    targetCellUnique =   xydat$occCellids [!xydat$occCellids %in% uniqueCellBuff]
    if (terra::nlyr(r.env)>1) df.CellTarget = terra::extract(r.env,targetCellUnique)
    if (terra::nlyr(r.env)==1) {
      df.CellTarget = terra::extract(r.env,targetCellUnique)
      df.CellTarget = data.frame (ID=1:length(df.CellTarget),LAYER = df.CellTarget)
      names (df.CellTarget) = c('ID',names(r.env))
    } 
    df.CellTarget$cellID = targetCellUnique
    df.cells = rbind (df.CellBuff,df.CellTarget)
    
    env.test = mymclapply (1:nrow(xydat), function (i){

      
      if (all(is.na(cellIds.directions[[i]]))) (return(NA))
      cellIds.directions[[i]] = stats::na.omit (cellIds.directions[[i]])
      NcellReps = table (cellIds.directions[[i]])
      
      df.env.cellBuffs  = lapply (1:length(NcellReps), function (m){
        cellID= as.numeric (names (NcellReps[m]))
        times = as.numeric (NcellReps[m])
        df.env.var = df.cells[which(df.cells$cellID ==cellID ), ]
        do.call("rbind", replicate(times, df.env.var, simplify = FALSE))
        
      })
      df.env.cellBuffs  = do.call(rbind,df.env.cellBuffs)
      df.env.cellBuffs  = df.env.cellBuffs[stats::complete.cases(df.env.cellBuffs),]
      df.env.targetcell = df.cells[which(df.cells$cellID == xydat$occCellids[i] ), ]
      
      if (nrow (df.env.cellBuffs)==0) {
        vars.outlier = rep(NA,length(names (r.env))+1 )
        names (vars.outlier) = c(names (r.env),'all')
        vars.outlier = as.data.frame (t(vars.outlier))
        return (vars.outlier)
      }
      
      if (any ( is.na (df.env.targetcell [,names (r.env)]) )) {
        vars.outlier = rep(NA,length(names (r.env))+1 )
        names (vars.outlier) = c(names (r.env),'all')
        vars.outlier = as.data.frame (t(vars.outlier))
        return (vars.outlier)
      }
      
      
      vars.outlier  <- sapply (names (r.env), function (nvar){
        
        qprobs = quantile (df.env.cellBuffs[,nvar],probs=env.quantiles)
        if (df.env.targetcell[,nvar] < qprobs[1]) {(return(1))}
        if (df.env.targetcell[,nvar] > qprobs[2]) {(return(1))}
        return (0)
      })
      vars.outlier =  as.data.frame (t(vars.outlier))
      vars.outlier$all <- mean (as.numeric(vars.outlier[1,]),na.rm=TRUE)
      vars.outlier
      
    })
    
    env.test = do.call (rbind, env.test)
    out$geoenvLowAccuracy_envDiff_value = env.test$all
    out$geoenvLowAccuracy_envDiff_test = (env.test$all > accept.threshold.env)
    out$geoenvLowAccuracy_envDiff_comments = paste('envQuantiles=',paste(env.quantiles,collapse = ','),';PercOutVariablesTh=',accept.threshold.env)
    
  }
  
  #write final score
  out$geoenvLowAccuracy_score <- .gimme.score (out)
  return (out)
}


# timeAccuracy ====
#' @title Flag occurrences with temporal inaccuracies
#' @description Detect records with low temporal accuracy
#' @param df data.frame of species occurrences
#' @param tf character. column name in df containing the dataset with the date/time where the species is recorded
#' @param iniTime character or Date. indicates the initial time point considered in year-month-day format. Defaults to NA.
#' @param endTime character or Date. indicates the time point considered in year-month-day format. Defaults to NA.
#' @param do logical. Should tests be performed? Default TRUE
#' @param method character. Vector of methods to be used. See details. Default 'all'
#' @param verbose logical. Print messages? Default FALSE
#' @return data.frame
#' @keywords Analysis 
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @details timeAccuracy function will implement different methods to assess the temporal scale of the parameters.\cr
#' @family analysis
#' @details timeAccuracy methods .\cr
#' Current implmented methods are:
#' 'noDate' : tests for an existing value in the date field. \cr
#' 'noDateFormatKnown' : tests for an known date format. \cr
#' 'outDateRange' : tests for dates within the iniTime and endTime interval.
#' @examples 
#' #see examples in vignetteXtra-occTest
#' #library devtools()
#' #install_github ('vignetteXtra-occTest')
#' @export
timeAccuracy  <- function (df,tf,iniTime=NA,endTime=NA,do=T,
                           method='all',
                           verbose=F)
  {
  #output results
  out = data.frame (
                    timeAccuracy_noDate_value =NA,
                    timeAccuracy_noDate_test =NA,
                    timeAccuracy_noDate_comments =NA,
                    
                    timeAccuracy_noDateFormatKnown_value =NA,
                    timeAccuracy_noDateFormatKnown_test =NA,
                    timeAccuracy_noDateFormatKnown_comments =NA,
                    
                    timeAccuracy_outDateRange_value =NA,
                    timeAccuracy_outDateRange_test =NA,
                    timeAccuracy_outDateRange_comments =NA,
                    
                    timeAccuracy_score=NA)[1:nrow (df),]
  
  row.names(x = out) <- NULL
  if (!do) {return (out)}

  #has time stamp? 
  if (any(method %in% c('noDate','all'))) {
    if (!is.null(tf)){
      vecTime = df[,tf]
      dfEmptyVals = data.frame (na=is.na(vecTime),
                                #null=is.null(as.character(vecTime)==NULL),
                                isEmpty = (as.character(vecTime)==''),
                                isEmptySpace = (as.character(vecTime)== '  '))
      
      vecTimeLogi = apply(dfEmptyVals,MARGIN = 1,FUN = function (x){any (x==TRUE)})
      out$timeAccuracy_noDate_value = vecTimeLogi * 1
      out$timeAccuracy_noDate_test = vecTimeLogi 
      out$timeAccuracy_noDate_comments = 'checked for NA,NULL, blankSpace and tab' 
    }
  }
  
  #has format stamp? 
  if (any(method %in% c('noDateFormatKnown','all'))) {
    if (!is.null(tf)){
      vecTime = df[,tf]
      if (all(is.na(vecTime))) {if (verbose) warning('all time stamps are NA')}
      if (!all(is.na(vecTime))){
        vecTime_date_formated <- anytime::anydate(vecTime)
        vecTime_guessed <- !is.na (vecTime_date_formated) & (!is.na(vecTime) | (vecTime!='')| (vecTime!='\t'))
        vecTime_timeStamp  <- ifelse((is.na(vecTime) | (vecTime=='')| (vecTime=='\t')),yes = NA,no = T) 
        out$timeAccuracy_noDateFormatKnown_value = !vecTime_guessed *1 * vecTime_timeStamp
        out$timeAccuracy_noDateFormatKnown_test = (!vecTime_guessed) * vecTime_timeStamp
        out$timeAccuracy_noDateFormatKnown_comments =NA  
      }
    }
  }
  
  #is within time range ?
  if (any(method %in% c('outDateRange','all'))) {
    if ( !is.null(tf) & (!any(is.na(c(iniTime,endTime)))) & (!any(is.null(c(iniTime,endTime)))) ){
      vecTime = df[,tf]
      if (all(is.na(vecTime))) {if (verbose) warning('all time stamps are NA')}
      if (!all(is.na(vecTime))){
        vecTime_date_formated <- anytime::anydate(vecTime)
        vecTime_guessed <- !is.na (vecTime_date_formated) & (!is.na(vecTime) | (vecTime!='')| (vecTime!='\t'))
        vecTime_timeStamp  <- ifelse((is.na(vecTime) | (vecTime=='')| (vecTime=='\t')),yes = NA,no = T) 
        it = lubridate::as_date(iniTime)
        et = lubridate::as_date(endTime)
        target_interval = lubridate::interval(it,et)
        `%within%` <- lubridate::`%within%`
        in_interval = vecTime_date_formated %within% target_interval
        out$timeAccuracy_outDateRange_value = !in_interval * 1
        out$timeAccuracy_outDateRange_test = !in_interval
        out$timeAccuracy_outDateRange_comments = paste('interval:',as.character (target_interval))
      }
    }
  }

  #write final score
  out$timeAccuracy_score <- .gimme.score (out)
  return (out)
}
# landUseSelect ====
#' @title Selection of records within a specified land use 
#' @description Selection of records within a specified land use 
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param method character. Select 'in' (default) when records are selected when inside specified land use types (type 'out' otherwise)
#' @param .points.crs crs character. Indicate coordinate reference system
#' @param ras.landUse raster. Land use raster with integer codes for land use classes.
#' @param .landUseCodes numeric. Vector of specified land use codes to either select records (method 'in') or to exclude records (method 'out')
#' @param do logical. Should range analysis be performed?
#' @param verbose logical. Print messages? Default  FALSE
#' @param output.dir character. Output directory
#' @return data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @export
#' @import sf

landUseSelect <- function (df,
                           xf,
                           yf,
                           method='in', #out would be the other one
                           .points.crs,
                           ras.landUse,
                           .landUseCodes,
                           do=TRUE, verbose=FALSE,output.dir){
  
  out <- data.frame (landUse_wrongLU_value=NA,
                     landUse_wrongLU_test=NA,
                     landUse_wrongLU_comments= NA,
                     landUse_score = NA)[1:nrow (df),]
  row.names(out) <- NULL
  
  if (!do) {return (out)}
  
  #start human influence analysis
  xydat <- df[,c(xf,yf)]
  
  #get species presence
  data.sp.pres <- as.data.frame (xydat)
  data.sp.pres <- sf::st_as_sf(data.sp.pres,coords=c(xf,yf),crs = .points.crs)
  
  #start with land use raster
  if(! inherits(ras.landUse,'RasterLayer')) {warning ('no raster of land use provided. Test not performed'); return (out)}
  
  #accommodate projections
  if (is.null (.points.crs)){sf::st_set_crs(data.sp.pres, value=terra::crs(ras.landUse)); warning ('ASSUMING Points and HumanRaster data have the same projection')}
  if (!is.null (.points.crs)) {sf::st_set_crs(data.sp.pres,value = .points.crs)}
  if (sf::st_crs(data.sp.pres) != terra::crs(ras.landUse)){
    if (is.na (terra::crs (ras.landUse)) ) {sf::st_set_crs (data.sp.pres,value = NA)  ; warning ('ASSUMING Points and HumanRaster data have the same projection')}
    if (!is.na (terra::crs (ras.landUse)) ) {data.sp.pres <- sf::st_transform(data.sp.pres, terra::crs(ras.landUse) )}
  }
  
  lu.sp.pres <-terra::extract(ras.landUse, y = data.sp.pres, cells=FALSE, xy=TRUE)
  row.id.lu.NA <- which(is.na(lu.sp.pres[,2]))
  if (length (row.id.lu.NA) != 0) {
    output.comments <- rep (NA,length =nrow (xydat))
    output.comments [row.id.lu.NA] <- paste0('No land cover available for this record;' )
  }
  
  lu.sp.value = lu.sp.pres[,2]
  
  if (method %in% c('in','all'))   {lu.sp.pres <- ( lu.sp.value %in% .landUseCodes)}
  if (method %in% c('out'))        {lu.sp.pres <- (!lu.sp.value %in% .landUseCodes)}
  
  out$landUse_wrongLU_value= !lu.sp.pres * 1
  out$landUse_wrongLU_test= !lu.sp.pres
  out$landUse_wrongLU_comments= paste(method,.landUseCodes)
  
  #compute score
  out$landUse_score <-  .gimme.score (out)
  
  
  return (out)
}





