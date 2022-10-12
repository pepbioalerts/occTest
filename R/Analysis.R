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
#' k <- data.frame (x=c(runif (n = 100),NA),y=c(runif (n = 100),1000))
#' filterMissing(k,xf='x',yf='y')
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
#' @descriptions checks for duplicated coordinates in the occurrence data frame
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
#' k <- data.frame (x=c(runif (n = 100),1000),y=c(runif (n = 100),1000),Reason=NA)
#' duplicatesexcludeAnalysis(k,xf='x',yf='y',resolution.in.minutes=60)
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
    ext <- raster::extent (raster.grid)
    resolution.raster <- raster::res (raster.grid)[1]
    rst <- raster::raster(xmn = ext@xmin, xmx = ext@xmax, ymn = ext@ymin, ymx = ext@ymax,
                  res = resolution.raster)
    }
  if(is.null(raster.grid)) {
    rst <- raster::raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90,
                  res = resolution.in.minutes/60)
  }

  #get dups in gridcell
  rst[] <- 1:raster::ncell(rst)
  spp <- as.character(unique(df$Species))
  xy <- data.frame(biogeo::coord2numeric(df[,xf]),
                   biogeo::coord2numeric(df[,yf]))
  cid <- raster::cellFromXY(rst, xy)
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
#' @descriptions Reassign coastal coordinates as needed to be in 
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
  ce0 <- raster::cellFromXY(rst, dd)
  vals <- raster::extract(rst, dd)
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
    ce1 <- raster::cellFromXY(rst, dd[f, ])
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
      a <- raster::adjacent(rst, ce3[i], directions = 8, pairs = FALSE,
                    target = NULL, sorted = FALSE, include = FALSE, id = FALSE)
      idx <- id2[i]
      # idx <- ce3[i] 
      b <- data.frame(i, a, id2 = idx)
      bb <- rbind(bb, b)
    }
    xy <- raster::xyFromCell(rst, bb$a, spatial = FALSE)
    vals <- raster::extract(rst, xy)
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
        dst <- raster::pointDistance(c(pts$x1, pts$y1), nce2, longlat = FALSE)
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
#' @param .points.proj4string Proj4string for the occurrence data
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
                                    .points.proj4string,
                                    .countries.shapefile,
                                    cfsf,
                                    excludeUnknownRanges = FALSE,
                                    excludeNotmatchCountry = FALSE,
                                    doRangeAnalysis=TRUE,
                                    verbose=FALSE) {
  

  if (!doRangeAnalysis) {
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
      if (is.null(.points.proj4string)) {
        .points.proj4string <- .countries.shapefile@proj4string
      if(verbose) print (paste ('ASSUMING points in projection',.countries.shapefile@proj4string))
      }
      sp.xydat <- sp::SpatialPoints(xydat,proj4string = .points.proj4string)
      overlay.sp.xydat <- sp::over(sp.xydat, .countries.shapefile)
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
#' @details
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param .c.field character. Country field in the species data.frame (df)
#' @param .ntv.ctry character. ISO3 country codes where species are considered native
#' @param idf character. Column with the taxon observation ID.
#' @param .inv.ctry character. ISO3 country codes where species are considered alien
#' @param .points.proj4string Proj4string argument. Coordinate reference system
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
                                .points.proj4string,
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
    #load centroid data
    dest_url = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/centroids.rds'
    outFile = paste0(tempdir(),'/centroids.rds')
    if (!file.exists(outFile)) utils::download.file(url=dest_url,destfile = outFile)
    centroid_data =  readRDS (outFile)
    
    #reformat dataframe of occurrences to be the same as in centroid Maitner methods
    if (is.null(idf)) { toid <- 1:nrow (df) } else {toid <- df[,idf]}
    occurrences.df = data.frame (taxonobservation_id=toid,
                                 country ='',
                                 state_province = '',
                                 county = '',
                                 latitude = df[,yf],
                                 longitude = df[,xf])
    
    #library(GNRS)
    #cleaned_countries<-GNRS_super_simple(country = df[,cf])
    
    #if (!is.null(cf)) {occurrences.df[,'country'] =  as.character(df[,cf])}
    #if ('state_province' %in% names (df)) {occurrences.df$state_province <- df$state_province}
    #if ('county' %in% names (df)) {occurrences.df$county <- df$county}
    #occurrences.df$cleaned_country<-NA
    
    #get admin unites from coordinates
    #### TO REPLACE BY A FUNCTION  .coords2country
    xydat <- df[,c(xf,yf)]
    if (!is.null(.countries.shapefile)){
      
      if (is.null(.points.proj4string)) {.points.proj4string <- .countries.shapefile@proj4string; if(verbose) print (paste ('ASSUMING points in projection',.countries.shapefile@proj4string))}
      sp.xydat <- sp::SpatialPoints(xydat,proj4string = .points.proj4string)
      overlay.sp.xydat <- sp::over(sp.xydat, .countries.shapefile)
      country_ext <- overlay.sp.xydat[,cfsf]
      occurrences.df$country <- country_ext
      
    }
    
    #extract countries of the points with when a country shapefile is not provided
    if (is.null(.countries.shapefile)){
      country_ext <-  .coords2country (xydat)
      occurrences.df$country <- country_ext
    }
    cleaned_countries<-try ({GNRS::GNRS_super_simple(country = country_ext)},silent=TRUE)
    if (inherits(cleaned_countries,'error') | is.null(cleaned_countries)){
      out$centroidDetection_BIEN_value <- NA
      out$centroidDetection_BIEN_test <- NA
      out$centroidDetection_BIEN_comments <- 'Not performed. GNRS error'
      
    } else{
      matchedCtry <- cleaned_countries$match_status=='full match'
      cleaned_countries <- cleaned_countries [matchedCtry,]
      maxValLoop = nrow(cleaned_countries)
      for(i in 1:nrow(cleaned_countries)){
        svMisc::progress(i,max.value = maxValLoop)
        occurrences.df$country <- as.character(occurrences.df$country)
        occurrences.df$country[which(occurrences.df$country ==cleaned_countries$country_verbatim[i])] <- unlist(cleaned_countries$country[i]) 
      }
      
      for(i in 1:nrow(cleaned_countries)){
        svMisc::progress(i,max.value = maxValLoop)
        occurrences.df$state_province <- as.character(occurrences.df$state_province)
        occurrences.df$state_province[which(occurrences.df$country ==cleaned_countries$country[i] &
                                              occurrences.df$state_province ==cleaned_countries$state_province_verbatim[i])] <- unlist(cleaned_countries$state_province[i]) 
      }
      
      for(i in 1:nrow(cleaned_countries)){
        svMisc::progress(i,max.value = maxValLoop)
        occurrences.df$county <- as.character(occurrences.df$county)
        occurrences.df$county[which(occurrences.df$country ==cleaned_countries$country[i] &
                                      occurrences.df$state_province ==cleaned_countries$state_province[i]&
                                      occurrences.df$county==cleaned_countries$county_parish_verbatim[i])] <- unlist(cleaned_countries$county_parish[i]) 
      }
      
      
      #Get relative distances
      k <- centroid_assessment(occurrences = occurrences.df, centroid_data = centroid_data)
      
      #convert relative distances to binary using threshold
      centroid_threshold <- 0.002 #smaller values require points to be closer to a centroid in order to be flagged
      
      k$is_centroid <- apply(X = k,MARGIN =  1,FUN = function(x){
        y<-as.vector(stats::na.omit(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative']))))
        #print(y)
        logiCentroid = any(as.numeric(y)<centroid_threshold)
        logiCentroid
      })
      k$is_centroid_val <- apply(X = k,MARGIN =  1,FUN = function(x){
        
        #print(x)
        if(all(is.na(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative'])))){
          
          return(NA) 
          
        }
        valOut = as.numeric(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative'])))[which.min(as.numeric(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative']))))]
        valOut
      })
      
      #construct comments
      k$comments <-apply(X = k,MARGIN = 1,FUN = function(x){
        #print(x)
        if(x['is_centroid']==1){ineq<-"<="}else{ineq<-">"}
        
        if(all(is.na(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative'])))){
          
          return("Centroid determination not conducted due to missing/unmatched political divison data" ) 
          
        }
        
        paste("Centroid determination based on ",
              c('country_cent_dist_relative','state_cent_dist_relative','county_cent_dist_relative')[which.min(as.numeric(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative']))))],
              ": ",
              as.numeric(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative'])))[which.min(as.numeric(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative']))))],
              ineq,
              "threshold value of",
              centroid_threshold
              
        )#paste
        
      } )#apply
      
      
      colnames(k)[which(colnames(k)=="is_centroid")]<-"centroidDetection_BIEN_test"
      colnames(k)[which(colnames(k)=="comments" )]<-"centroidDetection_BIEN_comments"
      out$centroidDetection_BIEN_value <- k$is_centroid_val
      out$centroidDetection_BIEN_test <- as.logical(k$centroidDetection_BIEN_test)
      out$centroidDetection_BIEN_comments <- k$centroidDetection_BIEN_comments
      
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
#' @param .points.proj4string character. Points coordinate projection.
#' @param ras.hii raster. Raster map of human influence index use
#' @param .th.human.influence numeric. Threhold to identify places of high human influence
#' @param .points.proj4string proj4string argument for df
#' @param do logical. Should range analysis be performed? Default TRUE
#' @param verbose logical. Print messages? Default FALSE
#' @param output.dir character. Output directory 
#' @return data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr). A Zizka (CoordinateCleaner functions)
#' @seealso \link[CoordinateCleaner]{cc_urb} for CoordinateCleaner functions
#' @family analysis
#' @export

humanDetection <- function (df,
                            xf,
                            yf,
                            method='all',
                            ras.hii,
                            .th.human.influence,
                            .points.proj4string,
                            do=TRUE, verbose=FALSE,output.dir){

  out <- data.frame (HumanDetection_HumanInfluence_value=NA,
                     HumanDetection_HumanInfluence_test=NA,
                     HumanDetection_HumanInfluence_comments= NA,
                     HumanDetection_UrbanAreas_value=NA,
                     HumanDetection_UrbanAreas_test=NA,
                     HumanDetection_UrbanAreas_comments=NA,
                     HumanDetection_score=NA
                     )[1:nrow (df),]
  row.names(out) <- NULL

  if (!do) {return (out)}

  #start human influence analysis
  xydat <- df[,c(xf,yf)]

  #get species presence
  data.sp.pres <- as.data.frame (xydat)
  sp::coordinates (data.sp.pres) <- stats::as.formula (paste0('~',xf,'+',yf))

  #start human influence index test
  if (any (method %in% c('hii','all'))) {

    class(ras.hii)=='RasterLayer'

    #accomodate projections
    if (is.null (.points.proj4string)){sp::proj4string(data.sp.pres) <- raster::projection(ras.hii); warning ('ASSUMING Points and HumanRaster data have the same projection')}
    if (!is.null (.points.proj4string)) {sp::proj4string(data.sp.pres) <- .points.proj4string}
    if (sp::proj4string(data.sp.pres) != raster::projection (ras.hii)){
      if (is.na (raster::projection (ras.hii)) ) {sp::proj4string(data.sp.pres) <- NA ; warning ('ASSUMING Points and HumanRaster data have the same projection')}
      if (!is.na (raster::projection (ras.hii)) ) {data.sp.pres <- sp::spTransform(data.sp.pres,CRSobj = raster::projection(ras.hii) )}
    }

    hii.sp.pres <-raster::extract(ras.hii, y = data.sp.pres, cellnumbers=FALSE, df=TRUE)
    row.id.hii.NA <- which(is.na(hii.sp.pres[,2]))
    if (length (row.id.hii.NA) != 0) {
      output.comments <- rep (NA,length =nrow (xydat))
      output.comments [row.id.hii.NA] <- paste0('No human filter available for this record;' )
    }
    
    hii.sp.value = hii.sp.pres[,2]

    hii.sp.pres <-( hii.sp.value >= .th.human.influence)

    out$HumanDetection_HumanInfluence_value=hii.sp.value
    out$HumanDetection_HumanInfluence_test=hii.sp.pres
    out$HumanDetection_HumanInfluence_comments= paste('Threshold influence=',.th.human.influence)

  }

  #start urban areas
  if (any (method %in% c('urban','all')))  {

    #check ref exists
    newoutdir = paste0(tempdir(),'/spatialData')
    alreadyDownloaded = file.exists (x = paste0(newoutdir,'/NE_urbanareas.shp'))
    if (alreadyDownloaded) myRef = rgdal::readOGR(dsn = newoutdir,layer = 'NE_urbanareas',verbose = FALSE)
    if (!alreadyDownloaded) myRef= NULL 
    
    cc_urb_test =  .cc_urb_occTest(x = df, lon =xf,lat=yf ,value='flagged', verbose = FALSE,ref = myRef,outdir = output.dir )
    cc_urb_test <- (!  cc_urb_test) * 1
    out$HumanDetection_UrbanAreas_value <- cc_urb_test
    out$HumanDetection_UrbanAreas_test <- as.logical(cc_urb_test) 
    out$HumanDetection_UrbanAreas_comments <- c('Urban areas from rnaturalearth')
  }
  
  #compute score
  out$HumanDetection_score <-  .gimme.score (out)


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
#' @param .points.proj4string proj4string argument for df
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
                                 .points.proj4string,
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
#' @descripion Detect geographical outliers using several tests
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param .alpha.parameter parameter setting for alphahull
#' @param .distance.parameter numeric. Maximum distance allowed. Default to 1000
#' @param .medianDeviation.parameter  numeric. Deviation parameter to . Default to 0.1
#' @param method charcter. Vector of methods to be used. See details. Default 'all'
#' @param .projString proj4string character. Indicate coordinate reference system
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
                                method = 'all',
                                .projString ,
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
                    geoOutliers_score=NA
                    )[1:nrow (df),]

  row.names (out) <- NULL

  if (!do) {return (out)}


  #prepare df for different methods of outliers in different packages
  df$species <- 'MyFakeSp'
  rownames(df) <-1:nrow(df)
  xydat <- df[,c(xf,yf)]
  if (nrow (xydat)<5) {warning('geoOutliers not perfomred. N<5'); return(out)}

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
                                            mltpl =  .medianDeviation.parameter ,
                                            sampling_thresh = .samplingIntensThreshold.parameter ,
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
    spdf = sp::SpatialPoints(coords = xydat,proj4string =.projString )
    if (length (spdf)<=5) {
      out$geoOutliers_Grubbs_value = rep(NA,times=nrow(out))
      out$geoOutliers_Grubbs_test = rep(NA,times=nrow(out))
      out$geoOutliers_Grubbs_comments = rep(paste('N<5. Analysis not run'),times=nrow(out))
      }
    
    if (length (spdf)>5){
      grubbs_tst = rep(0,times=nrow(out))
      grubbs.outl = try(occTest::findSpatialOutliers(myPres = spdf,verbose = FALSE),silent=TRUE)
      
      if (inherits(grubbs.outl,what = 'try-error')) {
        grubbs_tst = rep(NA,times=nrow(out)) 
        grubbs_comment <- rep (paste('Error in findSpatialOutliers'),times=nrow(out))
      }
      
      if (!inherits(grubbs.outl,what = 'try-error')) {
        grubbs_tst [grubbs.outl] = 1
        grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))
      }
      
      out$geoOutliers_Grubbs_value = grubbs_tst
      out$geoOutliers_Grubbs_test = as.logical(grubbs_tst)
      out$geoOutliers_Grubbs_comments = grubbs_comment
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
#' @param .projString proj4string character. Indicate coordinate reference system
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
#' @export
envOutliers  <- function (
                          df,
                          xf,
                          yf,
                          .r.env,
                          .th.perc.outenv ,
                          .sp.name ,
                          method='all',
                          .projString ,
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
                    envOutliers_score=NA
                    )[1:nrow (df),]

  if (!do) {return (out)}

  #prepare df for different methods of outliers in different packages
  df$species <- 'MyFakeSp'
  rownames(df) <-1:nrow(df)
  xydat <- df[,c(xf,yf)]
  
  if (nrow (xydat)<5) {
    warning('envOutliers not perfomred. N<5')
    
    return(out)}

  dat.environment <- raster::extract(x = .r.env, y= df[,c(xf,yf)], df = TRUE)
  
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
    environmental.columns <- 2:ncol(dat.environment)
    a.df <- lapply (environmental.columns, function (i){

      ev<-dat.environment[,i]

      try (a<-biogeo::outliers(rid=1:nrow(dat.environment), species, dups, ev), silent=TRUE)
      if (exists ('a') ) {
        a<- as.data.frame (a)
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}

      if (exists ('a') == FALSE ) {
        #in case there is no variation across the variable and then outliers function throws an error
        a <- data.frame (col1=rep (0,nrow(dat.environment)), col2 = rep (0,nrow(dat.environment) ))
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}



    })
    a.df <- do.call (cbind, a.df)
    choice.of.method.df <- a.df [,grep (names(a.df),pattern = outlier.method.biogeo)]

    if (inherits(.r.env,"RasterLayer")) {outlier.level <- choice.of.method.df}
    if (inherits(.r.env,"RasterStack")) {outlier.level <- round (rowMeans(choice.of.method.df) ,digits=2)}
    if (inherits(.r.env,"RasterBrick")) {outlier.level <- round (rowMeans(choice.of.method.df) ,digits=2)}

    outlier.env <- (outlier.level >= .th.perc.outenv)*1
    comments.outlier.env <- paste (paste("Outlier",outlier.method.biogeo,".level="), outlier.level)

    out$envOutliers_bxp_value = outlier.env
    out$envOutliers_bxp_test = as.logical(outlier.env)
    out$envOutliers_bxp_comments = comments.outlier.env

  }

  if (any (method %in% c('grubbs','all')))  {

    #spdf = sp::SpatialPointsDataFrame(data = dat.environment[,-1],coords = df[,c(xf,yf)] ,proj4string =.projString )

    env.grubbs_tst = rep(0,times=nrow(out))
    env.grubbs.outl = occTest::findEnvOutliers(myPres = dat.environment[,-1],myEnv = NULL,verbose = FALSE)
    env.grubbs_tst [env.grubbs.outl] = 1
    env.grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))
    
    
    out$envOutliers_Grubbs_value = env.grubbs_tst
    out$envOutliers_Grubbs_test = as.logical(env.grubbs_tst)
    out$envOutliers_Grubbs_comments = env.grubbs_comment

  }


  out$envOutliers_score <-  .gimme.score (out)
  return (out)
}

# geoEnvAccuracy ====
#' @title Coordinate accuracy
#' @description Detect records with low accuracy in space and time 
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param af character. column name in df containing the coordinate uncertainty value (in the same)
#' @param dsf character. column name in df containing the dataset to which the record belongs to (e.g. Forest Inventory of Spain)
#' @param ef character. column name in df containing the registered elevation for the record. 
#' @param tf character. column name in df containing the dataset with the date/time where the species is recorded
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
#' @export

geoEnvAccuracy  <- function (df,
                             xf,
                             yf,
                             af,dsf,ef,tf,
                             
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
                    
                    geoenvLowAccuracy_noDate_value =NA,
                    geoenvLowAccuracy_noDate_test =NA,
                    geoenvLowAccuracy_noDate_comments =NA,
                    
                    geoenvLowAccuracy_noDateFormatKnown_value =NA,
                    geoenvLowAccuracy_noDateFormatKnown_test =NA,
                    geoenvLowAccuracy_noDateFormatKnown_comments =NA,
                    
                    geoenvLowAccuracy_outDateRange_value =NA,
                    geoenvLowAccuracy_outDateRange_test =NA,
                    geoenvLowAccuracy_outDateRange_comments =NA,
                    
                    
                    
                    geoenvLowAccuracy_score=NA
  )[1:nrow (df),]
  
  row.names(x = out) <- NULL
  if (!do) {return (out)}
  #check if need parallel
  os =  .get_os()
  # if (doParallel==TRUE) {mymclapply <-  .hijack (parallel::mclapply,mc.cores=mc.cores)}
  # if (doParallel==TRUE & grepl(pattern = 'win',x = os)) {mymclapply <-  .hijack (parallelsugar::mclapply,mc.cores=mc.cores)}
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
        raster.elevation = try  (expr = {.getSRTM (xydat=df,download=TRUE, verbose=FALSE)},silent = TRUE)
        
      }
      #determining elevation test
      if (inherits(raster.elevation,"RasterLayer") ) {
        elev.xy <- raster::extract (raster.elevation, df[,c(xf,yf)])
        elev.difference <- abs (df[,ef] - elev.xy)
        out$geoenvLowAccuracy_elevDiff_value <- elev.difference
        out$geoenvLowAccuracy_elevDiff_test <- as.logical (ifelse (elev.difference >= elev.threshold, 1, 0))
        out$geoenvLowAccuracy_elevDiff_comments <- paste ('Elev difference:', elev.difference)
      }
      #skipping elevation test
      if (!inherits(raster.elevation,"RasterLayer")   ) {
        if(verbose) warning ("raster.elevation is not a raster; skipping elevDiff method")
      }
    }
    
  }
  
  #has time stamp? 
  if (any(method %in% c('noDate','all'))) {
    if (!is.null(tf)){
      vecTime = df[,tf]
      dfEmptyVals = data.frame (na=is.na(vecTime),
                                #null=is.null(as.character(vecTime)==NULL),
                                isEmpty = (as.character(vecTime)==''),
                                isEmptySpace = (as.character(vecTime)== '  '))
      
      vecTimeLogi = apply(dfEmptyVals,MARGIN = 1,FUN = function (x){any (x==TRUE)})
      out$geoenvLowAccuracy_noDate_value = vecTimeLogi * 1
      out$geoenvLowAccuracy_noDate_test = vecTimeLogi 
      out$geoenvLowAccuracy_noDate_comments = 'checked for NA,NULL, blankSpace and tab' 
    }
    
  }
  
  #has format stamp? 
  if (any(method %in% c('noDateFormatKnown','all'))) {
    if (!is.null(tf)){
      
      warning ('method noDateFormatKnown under development. Sorry')
      # vecTime = df[,tf]
      # 
      # dfVectime = data.frame(vecTime[1])
      # 
      # a = identify_dates(dfVectime)
      # dfVectime = data.table::as.data.table(dfVectime)
      # k = find_and_transform_dates(dfVectime)
      # 
      # 
      # out$geoenvLowAccuracy_noDate_value = vecTimeLogi * 1
      # out$geoenvLowAccuracy_noDate_test = vecTimeLogi 
      # out$geoenvLowAccuracy_noDate_comments = 'checked for NA,NULL, blankSpace and tab' 
    }
    
  }
  
  #is within time range ?
  if (any(method %in% c('outDateRange','all'))) {
    if (!is.null(tf)){
      
      warning ('method outDateRange under development. Sorry')

    }
    
  }
  
  #Need coordinate uncertainty analysis ? If not gimme score and leave
  if ( ! any (method %in% c('percDiffCell','envDeviation','all')) ) {
    #write final score
    out$geoenvLowAccuracy_score <- .gimme.score (out)
    return (out)}
  
  #Need coordinate uncertainty analysis ? Can you do it?
  if (is.null(af) | all(is.na(df[,af])) ) {
    print ('no coordinate accuracy/uncertainty field provided')
    #write final score
    out$geoenvLowAccuracy_score <- .gimme.score (out)
    return (out)
  }
  
  #Start methods requiring coordinate uncertainty
  if(length(af)>1) {df$new_accuracy = pmax(df[,af[1]],df[,af[2]],na.rm=TRUE); af = 'new_accuracy'}
  xydat = df[,c(xf,yf,af)]
  xydat$occCellids = raster::cellFromXY(r.env,as.data.frame (df[,c(xf,yf)]))
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
        dest.point.cellID = raster::cellFromXY(object = r.env, dest.point)
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
    if (raster::nlayers(r.env)>1) df.CellBuff = raster::extract(r.env,uniqueCellBuff,df=TRUE)
    if (raster::nlayers(r.env)==1) {
      df.CellBuff = raster::extract(r.env,uniqueCellBuff)
      df.CellBuff = data.frame (ID=1:length(df.CellBuff),LAYER = df.CellBuff)
      names (df.CellBuff) = c('ID',names(r.env))
    } 
    df.CellBuff$cellID = uniqueCellBuff
    
    targetCellUnique =   xydat$occCellids [!xydat$occCellids %in% uniqueCellBuff]
    if (raster::nlayers(r.env)>1) df.CellTarget = raster::extract(r.env,targetCellUnique,df=TRUE)
    if (raster::nlayers(r.env)==1) {
      df.CellTarget = raster::extract(r.env,targetCellUnique)
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


# landUseSelect ====
#' @title Selection of records within a specified land use 
#' @descriptions Selection of records within a specified land use 
#' @param df data.frame of species occurrences
#' @param xf character. column name in df containing the x coordinates
#' @param yf character. column name in df containing the y coordinates
#' @param method character. Select 'in' (default) when records are selected when inside specified land use types (type 'out' otherwise)
#' @param .points.proj4string proj4string character. Indicate coordinate reference system
#' @param ras.landUse raster. Land use raster with integer codes for land use classes.
#' @param .landUseCodes numeric. Vector of specified land use codes to either select records (method 'in') or to exclude records (method 'out')
#' @param do logical. Should range analysis be performed?
#' @param verbose logical. Print messages? Default  FALSE
#' @param output.dir character. Output directory
#' @return data.frame
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @export

landUseSelect <- function (df,
                           xf,
                           yf,
                           method='in', #out would be the other one
                           .points.proj4string,
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
  sp::coordinates (data.sp.pres) <- stats::as.formula (paste0('~',xf,'+',yf))
  
  #start with land use raster
  if(! inherits(ras.landUse,'RasterLayer')) {warning ('no raster of land use provided. Test not performed'); return (out)}
  
  #accommodate projections
  if (is.null (.points.proj4string)){sp::proj4string(data.sp.pres) <- raster::projection(ras.landUse); warning ('ASSUMING Points and HumanRaster data have the same projection')}
  if (!is.null (.points.proj4string)) {sp::proj4string(data.sp.pres) <- .points.proj4string}
  if (sp::proj4string(data.sp.pres) != raster::projection (ras.landUse)){
    if (is.na (raster::projection (ras.landUse)) ) {sp::proj4string(data.sp.pres) <- NA ; warning ('ASSUMING Points and HumanRaster data have the same projection')}
    if (!is.na (raster::projection (ras.landUse)) ) {data.sp.pres <- sp::spTransform(data.sp.pres,CRSobj =raster::projection(ras.landUse) )}
  }
  
  lu.sp.pres <-raster::extract(ras.landUse, y = data.sp.pres, cellnumbers=FALSE, df=TRUE)
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


# centroidsBIEN ====
#' @title Centroid assessment using BIEN methods
#'
#' @descriptions Detect centroids in occurrences data frame using BIEN methods
#' @details
#' @param occurrences data.frame with occurrence data
#' @param centroid_data data.frame with centroid coordinates for different administrative entities. 
#' @return data.frame
#' @keywords internal
#' @author Brian Maitner

centroid_assessment<-function(occurrences,centroid_data){
  
  occurrences$county<-as.character(occurrences$county)
  occurrences$state_province<-as.character(occurrences$state_province)
  occurrences$country<-as.character(occurrences$country)
  #message("This function assumes that occurrences are in wgs84 format, and then converts to whatever projection is used in the centroid data.  Only one projection at a time may be used in the centroid data, but multiple centroids per political division are fine.")
  
  ### NOT ACTIVATED BECAUSE WE HAVE DONE THIS BEFORE
  #deal with nonsense points by dropping those lines
  # if(length(  which(occurrences$latitude>90 | occurrences$latitude< -90) )>0  ){
  #   occurrences<-occurrences[-which(occurrences$latitude>90 | occurrences$latitude< -90),]}
  # if(length(  which(occurrences$longitude>180 | occurrences$longitude< -180) )>0){
  #   occurrences<-occurrences[-which(occurrences$longitude>180 | occurrences$longitude< -180),]}
  # if(any(is.na(occurrences$latitude) | (is.na(occurrences$longitude)))){
  #   occurrences<-occurrences[-which(is.na(occurrences$latitude) | (is.na(occurrences$longitude))),]  
  #   
  # }
  
  temp_points<-sp::SpatialPoints(coords = occurrences[c("longitude","latitude")],proj4string = sp::CRS("+init=epsg:4326"))
  temp_points<-sp::spTransform(x = temp_points,CRSobj = sp::CRS(projargs = unique(as.character(centroid_data$projection))) )
  
  
  occurrences$latitude <- temp_points@coords[,'latitude']
  occurrences$longitude <- temp_points@coords[,'longitude']  
  rm(temp_points)
  
  #output<-NULL
  
  output<-as.data.frame(matrix(nrow = nrow(occurrences),ncol=13))
  
  colnames(output)<-c("taxonobservation_id","country_cent_dist","country_cent_dist_relative","country_cent_type","country_max_uncertainty",
                      "state_cent_dist","state_cent_dist_relative","state_cent_type","state_max_uncertainty",
                      "county_cent_dist","county_cent_dist_relative","county_cent_type","county_max_uncertainty")
  output$taxonobservation_id<-occurrences$taxonobservation_id
  
  
  
  for(i in 1:length(unique(occurrences$country))){
    
    country_i<-unique(occurrences$country)[i]
    
    centroid_i<-centroid_data[which(centroid_data$gnrsed_country==country_i & is.na(centroid_data$state) & is.na(centroid_data$county)),]
    #test<-centroid_data[which(centroid_data$gnrsed_country==country_i),]
    
    if(nrow(centroid_i)>0){
      
      
      centroid_pt_i<-sp::SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      
      
      country_data_i<-occurrences[which(occurrences$country==country_i),]
      
      new_data<-matrix(nrow = nrow(country_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      country_pts_i<-sp::SpatialPointsDataFrame(coords = country_data_i[c("longitude","latitude")],data = as.data.frame(country_data_i$taxonobservation_id),proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- sp::spDistsN1(pts = country_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(country_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-sp::spDistsN1(pts = country_pts_i,pt = centroid_pt_i[p]) 
        
      }
      
      
      
      mins<-apply(X = new_data, MARGIN = 1,FUN = which.min)
      
      
      output$country_cent_dist_relative[which(occurrences$country==country_i)]<-unlist(apply(X = new_data,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      output$country_max_uncertainty[which(occurrences$country==country_i)]<-unlist(lapply(X = mins,FUN = function(x){centroid_i$max_uncertainty[x]}))
      output$country_cent_type[which(occurrences$country==country_i)]<-as.character(unlist(lapply(X = mins,FUN = function(x){centroid_i$centroid_type[x]}))  )
      
      output$country_cent_dist[which(occurrences$country==country_i)]<-unlist(apply(X = new_data_absolute,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      
      rm(centroid_pt_i)
      
      
    }
    rm(centroid_i)
    
  }#country level
  rm(country_data_i,new_data,country_i,country_pts_i,i,p,new_data_absolute,mins)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #State level
  
  
  for(i in 1:nrow(unique(occurrences[,c("country","state_province")]))){
    
    state_i<-unique(occurrences[,c("country","state_province")])[i,]
    
    centroid_i<-centroid_data[which(centroid_data$gnrsed_country==state_i$country & centroid_data$gnrsed_state==state_i$state_province & is.na(centroid_data$county)),]
    
    
    if(nrow(centroid_i)>0){
      
      centroid_pt_i<-sp::SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      
      state_data_i<-occurrences[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province),]
      
      new_data<-matrix(nrow = nrow(state_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      state_pts_i<-sp::SpatialPointsDataFrame(coords = state_data_i[c("longitude","latitude")],data = as.data.frame(state_data_i$taxonobservation_id),proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- sp::spDistsN1(pts = state_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(state_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-sp::spDistsN1(pts = state_pts_i,pt = centroid_pt_i[p]) 
        
      }
      
      
      
      mins<-apply(X = new_data, MARGIN = 1,FUN = which.min)
      
      output$state_cent_dist_relative[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province)]<-unlist(apply(X = new_data,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      output$state_max_uncertainty[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province)]<-unlist(lapply(X = mins,FUN = function(x){centroid_i$max_uncertainty[x]}))
      output$state_cent_type[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province)]<-as.character(unlist(lapply(X = mins,FUN = function(x){centroid_i$centroid_type[x]}))  )
      
      output$state_cent_dist[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province)]<-unlist(apply(X = new_data_absolute,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      
      rm(centroid_pt_i)
      
    }#if there is a matching centroid
    
    rm(centroid_i)
    
    
  }#state level
  #rm(state_data_i,new_data,state_i,state_pts_i,i,mins,p,new_data_absolute)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #county level
  
  for(i in 1:nrow(unique(occurrences[,c("country","state_province","county")]))){
    
    county_i<-unique(occurrences[,c("country","state_province","county")])[i,]
    
    centroid_i<-centroid_data[which(centroid_data$gnrsed_country==county_i$country & centroid_data$gnrsed_state==county_i$state_province & centroid_data$gnrsed_county==county_i$county  ),]
    
    
    if(nrow(centroid_i)>0){
      
      
      centroid_pt_i<-sp::SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      county_data_i<-occurrences[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county),]
      
      new_data<-matrix(nrow = nrow(county_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      county_pts_i<-sp::SpatialPointsDataFrame(coords = county_data_i[c("longitude","latitude")],data = as.data.frame(county_data_i$taxonobservation_id),proj4string = sp::CRS(projargs = unique(as.character(centroid_i$projection))))
      
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- sp::spDistsN1(pts = county_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(county_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data_absolute)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-sp::spDistsN1(pts = county_pts_i,pt = centroid_pt_i[p]) 
        
      }
      
      mins<-apply(X = new_data, MARGIN = 1,FUN = which.min)
      
      
      
      #Populate output file
      
      output$county_cent_dist_relative[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county)]<-unlist(apply(X = new_data,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      output$county_max_uncertainty[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county)]<-unlist(lapply(X = mins,FUN = function(x){centroid_i$max_uncertainty[x]}))
      output$county_cent_type[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county)]<-as.character(unlist(lapply(X = mins,FUN = function(x){centroid_i$centroid_type[x]}))  )
      output$county_cent_dist[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county)]<-unlist(apply(X = new_data_absolute,MARGIN = 1,FUN = function(x){x[which.min(x)]}) )
      
      rm(centroid_pt_i)
      
      
      
      
    }  
    rm(centroid_i)  
    
  }#county level
  
  #rm(county_data_i,new_data,county_i,county_pts_i,i,mins,p,new_data_absolute)
  
  return(output)
  
}




