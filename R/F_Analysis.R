#### FUNCTIONS FOR ANALYZING DATA IN THE WORKFLOW

#' @title Check for missing coordinates
#'
#' @description checks for missing coordinates in the occurrence dataframe
#' @details
#'
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @return list with two dataframes: stay = coordinates missing and continue = occurrence that you can retain for further analysis
#' @keywords internal
#'
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

filterMissing <- function (df, xf = x.field, yf = y.field, verbose=F){

 df.out <- df [!complete.cases  (df[,c(xf,yf)]),]
 df.continue <- df [complete.cases  (df[,c(xf,yf)]),]
 return (list (stay=df.out,continue=df.continue))

}

#' @title Duplicated records
#' @descriptions checks for duplicated coordinates in the occurrence dataframe
#' @details (inspired by biogeo duplicatesexclude function)
#'
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param resolution.in.minutes the resolution of environmental data used, specified in minutes
#' @param raster.grid An optional raster grid
#' @return list of three components: Dups.Exact = exact duplicate records, Dups.Grid= duplicates within environmental gridcell, continue = dataframe with good records
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


duplicatesexcludeAnalysis <- function (df=dat, xf=x.field, yf=y.field,
                                        resolution.in.minutes=res.in.minutes,
                                        raster.grid=NULL, verbose=F){


  #exact duplicates
  exact.dups <- duplicated (x = df[,c(xf,yf)]) *1
  df$Exclude <- exact.dups
  df$Reason [which (exact.dups==1)] <-'Exact Duplicated coordinates'
  df.exact.dups <- df[which(df$Exclude==1),]



  #update df
  df <- df[which(df$Exclude==0),]


  #build the grid for the duplicates
  if(!is.null(raster.grid)){
    ext <- extent (raster.grid)
    resolution.raster <- res (raster.grid)[1]
    rst <- raster(xmn = ext@xmin, xmx = ext@xmax, ymn = ext@ymin, ymx = ext@ymax,
                  res = resolution.raster)}

  if(is.null(raster.grid)) {
    rst <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90,
                  res = resolution.in.minutes/60)
  }

  #get dups in gridcell
  rst[] <- 1:ncell(rst)
  spp <- as.character(unique(df$Species))
  xy <- data.frame(biogeo::coord2numeric(df[,xf]),
                   biogeo::coord2numeric(df[,yf]))
  cid <- cellFromXY(rst, xy)

  dups <- (duplicated(cid)) * 1
  f1 <- which (dups==1)
  df$Exclude <- dups
  df$Reason <- as.character(df$Reason)
  df$Reason[f1] <- "Duplicated--GridCell"



  #indicate duplicates by grid cell (these needs updateing for efficiency)
  df.grid.dups <- df[which(df$Exclude==1),]


  #update df
  df <- df[which(df$Exclude==0),]


  #return
  return (list (Dups.Exact=df.exact.dups, Dups.Grid= df.grid.dups, continue = df))

}



#' @title SeaLand reassignement
#' @descriptions Reassign coastal coordinates as needed
#' @details (function inspiered in nearestcell in biogeo bu modified)
#' @param dat Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param rst An optional raster grid
#' @return list
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

nearestcell3 <- function (dat,
                          rst,
                          xf=x.field,
                          yf=y.field,
                          verbose=F) {

  dat$Correction <- as.character(dat$Correction)
  fx <- which(dat$Exclude == 0)
  x1 <- biogeo::coord2numeric(dat[,xf][fx])
  y1 <- biogeo::coord2numeric(dat[,yf][fx])
  datid <- dat$ID[fx]
  dd <- data.frame(x1, y1)
  ce0 <- raster::cellFromXY(rst, dd)
  vals <- raster::extract(rst, dd)

  #if the input raster is a dataframe
  if (class(vals)=='matrix'){
    vals <- apply(vals,1,FUN = sum)
    f <- which(is.na(vals))

  } else {
    f <- which(is.na(vals))
  }




  if (length(f) == 0){ if(verbose) print ("There are no missing values") ; return (dat)}

  if (length (f) != 0){
    id  <- datid[f]
    ce1 <- cellFromXY(rst, dd[f, ])

    if (any (is.na(ce1))) {

      if (length(ce1) == 1) {if(verbose) print("Coordinates in sea are out of environmental extent"); return(dat)}
      if(verbose) print("Some coordinates out of raster extent")
    }

    ff  <- which(!is.na(ce1))
    ce3 <- ce1[ff]
    dd2 <- dd[ff, ]
    id2 <- id[ff]
    bb  <- {}


    for (i in 1:length(ce3)) {
      a <- raster::adjacent(rst, ce3[i], directions = 8, pairs = FALSE,
                    target = NULL, sorted = FALSE, include = FALSE, id = FALSE)
      idx <- id2[i]
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
      g1 <- na.omit(g)
      near <- {
      }

      for (j in 1:length(uid)) {
        uj <- uid[j]
        fx <- which(g1$i == uj)
        nce <- g1[fx, ]
        nce2 <- as.matrix(cbind(nce$x, nce$y))
        pts <- dd2[uj, ]
        dst <- raster::pointDistance(c(pts$x1, pts$y1), nce2, longlat = F)
        fm <- which.min(dst)
        nr <- nce[fm, ]
        near <- rbind(near, nr)
      }

      mod <- format(Sys.time(), "%d-%m-%Y %H:%M:%S")
      dat$Modified <- as.character(dat$Modified)
      for (i in 1:nrow(near)) {
        f <- which(dat$ID == near$id2[i])
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
      datx <- list(dat = dat, moved = moved)
      return(datx)
    }


    else {
      warning("There are no records close enough to the nearest land/sea cells"); return(dat)
    }
  }
}


#' @title Range analysis
#'
#' @description Range analysis of species records based on invsive and native country
#' @details
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param .ntv.ctry Native country in ISO3
#' @param .inv.ctry Invasive country in ISO3
#' @param .c.field  Country field in the species dataframe
#' @param .points.proj4string Proj4string for dataframe
#' @param .countries.shapefile spatialPolygonDataFrame indicating countries
#' @param cfsf field of the country spatialPolygonDataFrame indicaing the country in ISO3
#' @return list
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

countryStatusRangeAnalysis=function(df=dat,
                                    xf=x.field,
                                    yf=y.field,
                                    .ntv.ctry=ntv.ctry,
                                    .inv.ctry=inv.ctry,
                                    .c.field=c.field,
                                    .points.proj4string=points.proj4string,
                                    .countries.shapefile=countries.shapefile,
                                    cfsf = countryfield.shapefile,
                                    excludeUnknownRanges = F,
                                    excludeNotmatchCountry = F,
                                    doRangeAnalysis=T,
                                    verbose=F) {

  if (!doRangeAnalysis) {
    df$countryStatusRangeAnalysis_wrongNTV_test   <- NA
    df$countryStatusRangeAnalysis_wrongCTRY_test  <- NA
    df$countryStatusRangeAnalysis_wrongINV_test   <- NA
    df$countryStatusRangeAnalysis_score <-NA

    stay = df[-1*1:nrow(df),]
    out <- list (stay=stay,continue = df)
    return (out)

  }


  # df=dat; xf=x.field; yf=y.field; .ntv.ctry=ntv.ctry;
  # .inv.ctry=inv.ctry;.c.field=c.field;.points.proj4string=points.proj4string;
  # .countries.shapefile=countries.shapefile;cfsf = countryfield.shapefile;exclude.unknown.ranges = F;
  # exclude.wrong.country.recorded = F; doRangeAnalysis=T; verbose=F

  #initial user information printing
  if (is.null (.ntv.ctry)  ) {if(verbose) print('WARNING: Species with no associated country. We assume all locations are native range')}
  if (is.null (.inv.ctry)  ) {if(verbose) print ('INFO: No invasive country provided')}
  if (is.null (.c.field)   ) {if(verbose) print ('INFO: No info on table of country of registration')}

  xydat <- df[,c(xf,yf)]

  #check country of coordinates
    #extract countries of the points with a given shapefile
    if (!is.null(.countries.shapefile)){

      if (is.null(.points.proj4string)) {.points.proj4string <- .countries.shapefile@proj4string; if(verbose) print (paste ('ASSUMING points in projection',.countries.shapefile@proj4string))}
      sp.xydat <- sp:::SpatialPoints(xydat,proj4string = .points.proj4string)
      overlay.sp.xydat <- sp:::over(sp.xydat, .countries.shapefile)
      fieldname  <- match(cfsf, names(overlay.sp.xydat))
      country_ext <- overlay.sp.xydat[, fieldname]


    }

    #extract countries of the points with when a country shapefile is not provided
    if (is.null(.countries.shapefile)){
      country_ext <- occProfileR:::.coords2country (xydat)
    }

  #inform for whether reported and country fo coordinates differ differ
  if (is.null (.c.field)) {
    if(verbose) print ("No country reported in occurrence database")
    wrong.ctry.reported <- rep (NA,length(country_ext))

  }
  if (!is.null (.c.field)) {

    #check if ISO3 format in country reported
    ctry.reported = as.character (df[,.c.field])




    #match reported country with extracted country
    wrong.ctry.reported <- (! as.character (df[,.c.field]) %in% as.character(country_ext))  * 1
  }

  #inform for which extracted countries are not in the native range
  if(is.null (.ntv.ctry)) {
    if(verbose) print ('No info of native country. Assuming all records are native')
    wrong.ntv.ctry.xy <- rep (NA,length(country_ext))
  }
  if(!is.null(.ntv.ctry) ) {wrong.ntv.ctry.xy <- (! country_ext %in% .ntv.ctry) * 1}

  #inform for which extracted countries are not in the invasive range
  if (is.null (.inv.ctry))   {
    if(verbose) print ('No info of invasive countries')
    wrong.inv.ctry.xy <- rep (NA,length(country_ext))
  }
  if (!is.null (.inv.ctry))  {wrong.inv.ctry.xy <- (! country_ext %in% .inv.ctry) * 1}

  #add to df
  df$countryStatusRangeAnalysis_wrongNTV_test   <- wrong.ntv.ctry.xy
  df$countryStatusRangeAnalysis_wrongCTRY_test <- wrong.ctry.reported
  df$countryStatusRangeAnalysis_wrongINV_test   <- wrong.inv.ctry.xy

  #Exclusion from subsequent analysis for wrong native and invasive ranges
  #  (depending how much data has been provided)
  if (excludeUnknownRanges == T & !is.null(.ntv.ctry)  & !is.null(.inv.ctry) ){
    total.data <- (!is.na (wrong.ntv.ctry.xy) *1) + (!is.na(wrong.inv.ctry.xy)*1)
    sum.wrong.extracted.ranges <-  mapply (function (x,y) (sum (x,y,na.rm = T)),wrong.ntv.ctry.xy, wrong.inv.ctry.xy)
    exclusion.based.on.wrong <- sum.wrong.extracted.ranges / total.data
    df$Exclude  <- (exclusion.based.on.wrong == 1)*1
    exclude.rows <- which (df$Exclude == 1)
    #add the reason of exclusion
    if (length(exclude.rows) > 0) {df$Reason [exclude.rows] <-"XY in not in invasive or native ranges"}
  }
  if (excludeUnknownRanges == F) {df$Exclude  <- 0}

  # Exclusion based on wrong country reported
  if  (excludeNotmatchCountry == T & !is.null (.c.field)) {

    exclusion.based.on.wrong <- df$Exclude  + (10*wrong.ctry.reported)
    only.wrong.reported.range <- which (exclusion.based.on.wrong == 10 )
    wrong.reported.range.and.informedrange <- which (exclusion.based.on.wrong == 11 )

    df$Reason [only.wrong.reported.range]  <- "Incongruence reported country and coordinates"
    df$Exclude [only.wrong.reported.range] <- 1
    df$Reason [wrong.reported.range.and.informedrange] <-"XY in not in invasive or native ranges;Incongruence reported country and coordinates"

  }

  #output
  df$countryStatusRangeAnalysis_score = occProfileR:::.gimme.score (df)
  out <- list (stay = df[which(df$Exclude==1),], continue = df[which(df$Exclude!=1),])
  return (out)

}


#' @title Centroid detection function
#'
#' @descriptions Detect centroids in occurrences dataframe
#' @details
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param .ntv.ctry Native country in ISO3
#' @param .inv.ctry Invasive country in ISO3
#' @param .points.proj4string Proj4string argument
#' @param .r.env R env
#' @return list
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

centroidDetection <- function (df=dat,
                                xf=x.field,
                                yf=y.field,
                                cf =c.field,
                                idf = taxonobservation.id,
                                .ntv.ctry=ntv.ctry,
                                .inv.ctry=inv.ctry,
                                .points.proj4string=points.proj4string,
                                .r.env=r.env,
                                .countries.shapefile=countries.shapefile,
                                 method='all',
                                 do = T, verbose=F){

  # df=dat;xf=x.field;yf=y.field;cf =c.field;
  # .ntv.ctry=ntv.ctry; .inv.ctry=inv.ctry;
  # .points.proj4string=points.proj4string
  # .r.env=r.env;method='all';do = T



  #table of outputs
  out = data.frame (centroidDetection_BIEN_test=NA,
                    centroidDetection_BIEN_comments =NA,
                    centroidDetection_speciesGeoCodeR_test=NA,
                    centroidDetection_speciesGeoCodeR_comments=NA,
                    centroidDetection_CoordinateCleaner_test=NA,
                    centroidDetection_CoordinateCleaner_comment=NA,
                    centroidDetection_score=NA
                    )[1:nrow (df),]

  row.names(out) <- NULL

  if(!do) { return (out)}
  

  #Method BIEN
  if (any(method %in% c('BIEN','all'))){
    
    #load centroid data
    centroid_data = system.file('ext/Centroids/centroids_equal_area_km_box_shape.csv',package='occProfileR')
    centroid_data =  read.csv (centroid_data)
    
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
    xydat <- df[,c(xf,yf)]
    if (!is.null(.countries.shapefile)){
      
      if (is.null(.points.proj4string)) {.points.proj4string <- .countries.shapefile@proj4string; if(verbose) print (paste ('ASSUMING points in projection',.countries.shapefile@proj4string))}
      sp.xydat <- sp:::SpatialPoints(xydat,proj4string = .points.proj4string)
      overlay.sp.xydat <- sp:::over(sp.xydat, .countries.shapefile)
      country_ext <- overlay.sp.xydat$name
      occurrences.df$country <- country_ext
      
    }
    
    #extract countries of the points with when a country shapefile is not provided
    if (is.null(.countries.shapefile)){
      country_ext <- occProfileR:::.coords2country (xydat)
      occurrences.df$country <- country_ext
    }
    
    cleaned_countries<-GNRS:::GNRS_super_simple(country = country_ext)
    
    
    for(i in 1:nrow(cleaned_countries)){
      occurrences.df$country[which(occurrences.df$country ==cleaned_countries$country_verbatim[i])] <- unlist(cleaned_countries$country[i]) 
    }
    
    for(i in 1:nrow(cleaned_countries)){
      occurrences.df$state_province[which(occurrences.df$country ==cleaned_countries$country[i] &
                                            occurrences.df$state_province ==cleaned_countries$state_province_verbatim[i])] <- unlist(cleaned_countries$state_province[i]) 
    }
    
    for(i in 1:nrow(cleaned_countries)){
      occurrences.df$county[which(occurrences.df$country ==cleaned_countries$country[i] &
                                    occurrences.df$state_province ==cleaned_countries$state_province[i]&
                                    occurrences.df$county==cleaned_countries$county_parish_verbatim[i])] <- unlist(cleaned_countries$county_parish[i]) 
    }
    
    
    #Get relative distances
    k <- centroid_assessment(occurrences = occurrences.df, centroid_data = centroid_data)
    
    #convert relative distances to binary using threshold
    centroid_threshold <- 0.002 #smaller values require points to be closer to a centroid in order to be flagged
    
    k$is_centroid <- apply(X = k,MARGIN =  1,FUN = function(x){
      y<-as.vector(na.omit(unlist(c(x['country_cent_dist_relative'],x['state_cent_dist_relative'],x['county_cent_dist_relative']))))
      #print(y)
      if(any(as.numeric(y)<centroid_threshold)){return(1)}else(return(0))
      
    })
    
    
    #construct comments
    k$comments <-
      
      apply(X = k,MARGIN = 1,FUN = function(x){
        #print(x)
        if(x['is_centroid']==1){ineq<-"<"}else{ineq<-">"}
        
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
        
      }
      )#apply
    
    
    colnames(k)[which(colnames(k)=="is_centroid")]<-"centroidDetection_BIEN_test"
    colnames(k)[which(colnames(k)=="comments" )]<-"centroidDetection_BIEN_comments"
    
    out$centroidDetection_BIEN_test <- k$centroidDetection_BIEN_test
    out$centroidDetection_BIEN_comments <- k$centroidDetection_BIEN_comments
    
    
  }
  #Method SpeciesGeoCodeR
  if (any(method %in% c('speciesGeoCodeR','all'))){
    if(!"speciesgeocodeR" %in% rownames(installed.packages())){
      message('You must first install the speciesGeoCodeR package. It does not yet work with R3.6 so you have to get it from https://github.com/azizka/speciesgeocodeR/' )
      return()
    }
    #load country reference
    #require (speciesgeocodeR)
    #data('countryref')
    countryref <- CoordinateCleaner::countryref
    #get target centroids countries and their captials
    country.names <- countryref$iso3

    #chunk if we were to match the countries provided by the user
    #country.names <- c(.ntv.ctry,.inv.ctry)
     #if (is.null(country.names)) {}

    centroid.ctry <- countryref [which (countryref$iso3 %in% country.names),c('centroid.lon','centroid.lat')]
    names (centroid.ctry) <- c('lon','lat')
    centroid.ctry.cap <-countryref [which (countryref$iso3 %in% country.names),c('capital.lon','capital.lat')]
    names (centroid.ctry.cap) <- c('lon','lat')
    centroids <- rbind (centroid.ctry,centroid.ctry.cap)
    centroids <- centroids [complete.cases(centroids),]

    #way out if no countires
    if (nrow (centroids) == 0 ) {
      if(verbose) print (paste ('No centroids found for', country.names))
      centroidDetection_speciesGeoCodeR_test <- 0
      centroidDetection_speciesGeoCodeR_comments <- 'No country nor capital centroid found'
      }

    #if we found them
    if  (nrow (centroids) >0 ) {
      #get raster cells in centroids
      cc <- centroids
      sp::coordinates (cc) <- ~lon+lat
      sp::proj4string(cc) <- .points.proj4string

      #deal with projections
      if (sp::proj4string(cc)  != raster::projection (.r.env)){
        if (is.na(projection(.r.env)) ) {
          sp::proj4string(cc) <- raster::projection(.r.env)
          warning ('ASSUMING Centroids and Environmental data have the
                   same projection')}
        if (!is.na(projection(.r.env)) ) {
          cc <- sp::spTransform(cc,CRSobj =projection(.r.env) )}
      }

      #get raster cells in centroids
      focal.cell.ids <- raster::extract (.r.env,y = cc, cellnumbers=T, df=T)
      focal.cell.ids <- focal.cell.ids$cells
      focal.cells.ids <- focal.cell.ids[is.na(focal.cell.ids)]
      neig.cells <- 8 #we could choose other
      adj.cells  <- raster::adjacent (.r.env,cells = focal.cell.ids, directions = neig.cells, include = F, id=F)
      cell.ids <- unique (as.vector(adj.cells))
      cell.ids <- cell.ids[!is.na(cell.ids)]
      if (length(cell.ids) == 0) cell.ids = 0

      #get species presence
      data.sp.pres <- as.data.frame (df[,c(xf,yf)])
      sp:::coordinates (data.sp.pres) <- as.formula (paste0('~',xf,'+',yf))

      #accomodate projections
      if (is.null (.points.proj4string)){sp::proj4string(data.sp.pres) <- raster::projection(.r.env); warning ('ASSUMING Centroids and Environmental data have the same projection')}
      if (!is.null (.points.proj4string)) {sp::proj4string(data.sp.pres) <- .points.proj4string}
      if (sp::proj4string(data.sp.pres) != raster::projection (.r.env)){
        if (is.na(raster::projection(.r.env))) {sp::proj4string(data.sp.pres) <- NA ; warning ('ASSUMING Centroids and Environmental data have the same projection')}
        if (!is.na (raster::projection (.r.env)) ) {data.sp.pres <- sp::spTransform(data.sp.pres,CRSobj =raster::projection(.r.env) )}
      }

      #get species presence cells
      cell.ids.sp <- raster::extract (.r.env,y = data.sp.pres, cellnumbers=T, df=T)

      #species presence in centroids or neighbors?
      sp.in.centroid <- cell.ids.sp$cells %in% cell.ids * 1

      out$centroidDetection_speciesGeoCodeR_test <- sp.in.centroid
      out$centroidDetection_speciesGeoCodeR_comments <- NA
    }
  }
  #Method CoordinateCleaner
  if (any(method %in% c('CoordinateCleaner','all'))){

    cc_cap_test = CoordinateCleaner::cc_cap(x = df,lon =xf,lat=yf,value='flagged', verbose = F)
    cc_cap_test = (!cc_cap_test) *1
    cc_cen_test <- CoordinateCleaner::cc_cen (x = df, lon =xf,lat=yf ,value='flagged', verbose = F)
    cc_cen_test <- (!cc_cen_test) * 1
    out$centroidDetection_CoordinateCleaner_test = ifelse ((cc_cap_test==1 | cc_cen_test==1),1,0)

  }

  out$centroidDetection_score = occProfileR:::.gimme.score (x = out)
  return (out)

}


#' @title HYPER HUMAN ENVIRONMENT FUNCTION
#'
#' @description Detect occurrences in heavily human-impacted environments
#' @details
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param .points.proj4string proj4string argument for dataframe
#' @param ras.hii Raster of human influence index
#' @param .th.human.influence threshold of human influence index
#' @return list
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

hyperHumanDetection <- function (df=dat,
                                  xf=x.field,
                                  yf=y.field,
                                 method='all',
                                  .points.proj4string=points.proj4string,
                                  ras.hii=ras.hii,
                                  .th.human.influence =th.human.influence,
                                  do=T, verbose=F){

  out <- data.frame (HumanInfluence_test=NA,
                     HumanInfluence_comments= NA,
                     UrbanAreas_test=NA,
                     UrbanAreas_comments=NA,
                     hyperHumanDetection_score=NA
                     )[1:nrow (df),]
  row.names(out) <- NULL

  if (!do) {return (out)}

  #start human influence analysis
  xydat <- df[,c(xf,yf)]

  #get species presence
  data.sp.pres <- as.data.frame (xydat)
  coordinates (data.sp.pres) <- as.formula (paste0('~',xf,'+',yf))

  #start human influence index test
  if (any (method %in% c('hii','all'))) {

    class(ras.hii)=='RasterLayer'

    #accomodate projections
    if (is.null (.points.proj4string)){proj4string(data.sp.pres) <- projection(ras.hii); warning ('ASSUMING Points and HumanRaster data have the same projection')}
    if (!is.null (.points.proj4string)) {proj4string(data.sp.pres) <- .points.proj4string}
    if (proj4string(data.sp.pres) != projection (ras.hii)){
      if (is.na (projection (ras.hii)) ) {proj4string(data.sp.pres) <- NA ; warning ('ASSUMING Points and HumanRaster data have the same projection')}
      if (!is.na (projection (ras.hii)) ) {data.sp.pres <- spTransform(data.sp.pres,CRSobj =projection(ras.hii) )}
    }

    hii.sp.pres <-raster::extract(ras.hii, y = data.sp.pres, cellnumbers=F, df=T)
    row.id.hii.NA <- which(is.na(hii.sp.pres[,2]))
    if (length (row.id.hii.NA) != 0) {
      output.comments <- rep (NA,length =nrow (xydat))
      output.comments [row.id.hii.NA] <- paste0('No human filter available for this record;' )
    }

    hii.sp.pres <-( hii.sp.pres[,2]>=.th.human.influence)*1

    out$HumanInfluence_test=hii.sp.pres
    out$HumanInfluence_comments= paste('Threshold influence=',.th.human.influence)

  }



  #start uban areas

  if (any (method %in% c('urban','all')))  {

    cc_urb_test = CoordinateCleaner::cc_urb(x = df, lon =xf,lat=yf ,value='flagged', verbose = F )
    cc_urb_test <- (!  cc_urb_test) * 1
    out$UrbanAreas_test <- cc_urb_test
    out$UrbanAreas_comments <- c('Urban areas from rnaturalearth')

    #get output
    out$hyperHumanDetection_score <- occProfileR:::.gimme.score (out)


  }


  return (out)
}


#' @title POTENTIAL BOTANIC GARDEN LOCALITY FROM NAME
#'
#' @description Detect occurrences potentially in botanical gardens via locality name
#' @details
#' @param df Data.frame of species occurrences
#' @param lf The locality field in the dataframe, df.
#' @return list
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
institutionLocality <- function (df=dat,
                                 xf=x.field,
                                 yf=y.field,
                                 .points.proj4string=points.proj4string,
                                 lf=l.field,
                                 method='all',
                                 do=T, verbose=F
                                 ){

  #output table
  out = data.frame (institutionLocality_fromBotanicLocalityName_test=NA,
                    institutionLocality_fromBotanicLocalityName_comments=NA,
                    institutionLocality_fromCoordinates_test=NA,
                    institutionLocality_fromCoordinates_comments=NA,
                    institutionLocality_score=NA
                    )[1:nrow (df),]

  row.names (out) <- NULL

  if (!do) {return (out)}

  if (any (method %in% c('fromBotanicLocalityName','all'))) {

    #if no locality data is provided
    if(is.null (lf)) {
      if(verbose) print ('No locality in input data')
      return (out)
    }


    #if locality data is not provided
    loc=df[,lf]
    potential.names <- c('botanic', 'botanische', 'botanico', 'jardin', 'garden', 'botanical')
    localities <- as.character (loc)
    bot.garden <-sapply (localities, function (x) {

      if (is.na(x)){return (0)}
      if (x==''){return (0)}

      m <- tolower(.multiple.strsplit(x,multiple.splits = c(' ',',',';',':','/','\\\\','[.]')))
      potential.bot.garden <- sum ((m %in% potential.names) *1)
      potential.bot.garden <- ifelse(potential.bot.garden>0, 1, 0)
      return (potential.bot.garden)

    })
    bot.garden <- as.numeric (bot.garden)


    out$institutionLocality_fromBotanicLocalityName_test = bot.garden
    out$institutionLocality_fromBotanicLocalityName_comments = 'likely botanical garden'


  }

  if (any (method %in% c('fromCoordinates','all'))) {
    cc_inst_test = CoordinateCleaner::cc_inst(x = df, lon =xf,lat=yf ,value='flagged', verbose = F )
    cc_inst_test <- (!  cc_inst_test) * 1
    cc_gbif_test = CoordinateCleaner::cc_gbif(x = df, lon =xf,lat=yf ,value='flagged', verbose = F )
    cc_gbif_test <- (!  cc_gbif_test) * 1

    cc_gbif_comments = ifelse(cc_gbif_test==1,'GBIFheadquarters',NA)
    cc_inst_comments = ifelse(cc_inst_test==1,'BiodivInst',NA)

    out$institutionLocality_fromCoordinates_test = ifelse ((cc_gbif_test==1 | cc_inst_test==1),1,0)
    out$institutionLocality_fromCoordinates_comments = occProfileR:::.paste3(cc_gbif_comments,cc_inst_comments)

  }

  out$institutionLocality_score <- occProfileR:::.gimme.score (out)

  return (out)

}


#' @title ALPHA HULL OUTLIERS
#'
#' @descripion Detect geographical outliers using alphahulls
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param .alpha.parameter parameter setting for alphahull
#' @return list
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
geoOutliers         <- function (df=dat,
                                xf=x.field,
                                yf=y.field,
                                .alpha.parameter=alpha.parameter,
                                .distance.parameter=1000,
                                .medianDeviation.parameter=5,
                                .samplingIntensThreshold.parameter=0.1,
                                .projString = points.| envOutliers    | 'bxp'  |  proj4string,
                                method = 'all',
                                do=T, verbose=F){
  

  #build dataframe of all potential results
  out = data.frame (geoOutliers_alphaHull_test=NA,
                    geoOutliers_alphaHull_comments = NA,
                    geoOutliers_distance_test=NA,
                    geoOutliers_distance_comments = NA,
                    geoOutliers_median_test=NA,
                    geoOutliers_median_comments=NA,
                    geoOutliers_quantileSamplingCorr_test=NA,
                    geoOutliers_quantileSamplingCorr_comments=NA,
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

  if (any (method %in% c('alphaHull','all'))){




    #We typically will consider by default an alpha 2 -like ALA
    if (nrow (df) > 5  & !is.na(.alpha.parameter)) {
      ah.2alpha <- alphahull::ahull (x = xydat[,xf], y= xydat[,yf], alpha=.alpha.parameter)

      points.outside.alphahull <- !alphahull::inahull(ah.2alpha,as.matrix(xydat))
      points.outside.alphahull <- points.outside.alphahull * 1
      #ahshape <- occProfileR:::.ah2sp(ah.2alpha)
      out.comments <- paste0('GeoIndicator alphaHull (Alpha=',.alpha.parameter,')')
    }


    if (nrow (df) <=5 ) {
      points.outside.alphahull <- rep (NA, nrow(df) )
      out.comments <- paste0('Not enough samples for GeoIndicator alphaHull.')
    }

    #write results into the results dataframe
    out$geoOutliers_alphaHull_test <- points.outside.alphahull
    out$geoOutliers_alphaHull_comments <- out.comments

    #This is a bit experimental, when the user does not provide a predetermined value of alpha
    # we perform some analysis to get an "optimal" alpha parameter.
    # not implemented yet
    # if (nrow (df) >5 & is.na(.alpha.parameter) ) {
    #
    #   #build alpha hulls and determine alpha parameter for each from 1 to 10
    #   ah.hulls.list <- lapply (1:10, function (a){
    #
    #     try(ah <- alphahull::ahull (x = df[,xf], y= df[,yf], alpha=a), silent = T)
    #     if (exists("ah")){return (ah)}
    #
    #   })
    #
    #   ah.hulls.thresholds <- lapply (1:10, function (i) {
    #     if((is.null(ah.hulls.list[[i]]) == F) & (is.list (ah.hulls.list[[i]]) == T)  ) {inahull(ah.hulls.list[[i]] , as.matrix(xydat) ) * i}
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
    #       m <- min (x, na.rm=T)
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
    #       if (exists("ah") & is.null(ah) == F & is.list(ah)==T ){
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
    #       alpha.choice <- alpha.choice[order(alpha.choice$optim.alpha.values,decreasing = T,na.last = T),]
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
    cc_dist_test = occProfileR::Mycc_outl(x = df,lon = xf,lat = yf,method = 'distance',value = 'flagged',tdi = .distance.parameter, verbose=F)
    cc_dist_test <- (!  cc_dist_test) * 1
    out$geoOutliers_distance_test  = cc_dist_test
    out$geoOutliers_distance_comments = rep(paste('dist.threshold=',.distance.parameter,'m'),times=nrow(out))
    if (nrow (df)<7) {out$geoOutliers_distance_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))}
  }
  

  if (any (method %in% c('median','all'))){

    cc_med_test = occProfileR::Mycc_outl(x = df,lon = xf,lat = yf,method = 'mad',value = 'flagged',mltpl =  .medianDeviation.parameter , verbose=F)
    cc_med_test <- (!  cc_med_test) * 1
    out$geoOutliers_median_test  = cc_med_test
    out$geoOutliers_median_comments = rep(paste('medDeviation=', .medianDeviation.parameter),times=nrow(out))
    if (nrow (df)<7) {out$geoOutliers_median_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))}
    
  }

  if (any (method %in% c('quantSamplingCorrected','all'))){
    

    cc_qsc_tst = occProfileR::Mycc_outl(x = df,lon = xf,lat = yf,method='quantile',value = 'flagged',
                                            mltpl =  .medianDeviation.parameter ,
                                            sampling_thresh = .samplingIntensThreshold.parameter ,
                                            verbose=F)
    cc_qsc_tst <- (!  cc_qsc_tst) * 1


    out$geoOutliers_quantileSamplingCorr_test = cc_qsc_tst
    out$geoOutliers_quantileSamplingCorr_comments = rep(paste('IQRmultiplier=', .medianDeviation.parameter,';SamplIntensityThres=',.samplingIntensThreshold.parameter),times=nrow(out))
    if (nrow (df)<7) {out$geoOutliers_quantileSamplingCorr_comments = rep(paste('N<7. Analysis not run'),times=nrow(out))}
    

  }

  if (any (method %in% c('grubbs','all'))) {
    spdf = sp::SpatialPoints(coords = xydat,proj4string =.projString )

    grubbs_tst = rep(0,times=nrow(out))
    grubbs.outl = occProfileR:::findSpatialOutliers(myPres = spdf,verbose = F)
    grubbs_tst [grubbs.outl] = 1
    grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))

    out$geoOutliers_Grubbs_test = grubbs_tst
    out$geoOutliers_Grubbs_comments = grubbs_comment
  }

  out$geoOutliers_score <- occProfileR:::.gimme.score (out)

  return (out)
}


#Fut Develop: Potentially, we could use bicolim, my little playing around has found that it acutally takes some time to do the biclim algorithm rather than rev.jacknife
#' @title ENVIRORNMENTAL OUTLIERS (using both jacknife and boxplot)
#'
#' @description Detect environmental outliers using jacknife and boxplot
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param .r.env R enviroment
#' @param .th.perc.outenv Something
#' @param .sp.name Species name
#' @param outlier.method  Method to use for detecting outliers.  Defaults to 'all'
#' @return list
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
envOutliers  <- function (.r.env=r.env,
                          df= dat,
                          xf=x.field,
                          yf =y.field,
                          .th.perc.outenv = th.perc.outenv,
                          .sp.name = sp.name,
                          method='all',
                          .projString = points.proj4string,
                          do=T, verbose=F){

  #output results
  out = data.frame (envOutliers_missingEnv_test = NA,
                    envOutliers_bxp_test = NA,
                    envOutliers_bxp_comments = NA,
                    envOutliers_Grubbs_test=NA,
                    envOutliers_Grubbs_comments=NA,
                    envOutliers_score=NA
                    )[1:nrow (df),]

  if (!do) {return (out)}

  #prepare df for different methods of outliers in different packages
  df$species <- 'MyFakeSp'
  rownames(df) <-1:nrow(df)
  xydat <- df[,c(xf,yf)]

  dat.environment <- extract(x = .r.env, y= df[,c(xf,yf)], df = T)
  missingEnvironment <- (!complete.cases(dat.environment))*1
  out$envOutliers_missingEnv_score <- missingEnvironment

  if (any (method %in% c('bxp','all'))) {

    outlier.method.biogeo = '-bxp'
    dat.environment <- extract(x = .r.env, y= df[,c(xf,yf)], df = T)
    missingEnvironment <- (!complete.cases(dat.environment))*1
    species <- rep(.sp.name, nrow(dat.environment))
    dups <- rep(0, nrow(dat.environment))
    environmental.columns <- 2:ncol(dat.environment)
    a.df <- lapply (environmental.columns, function (i){

      ev<-dat.environment[,i]

      try (a<-biogeo::outliers(rid=1:nrow(dat.environment), species, dups, ev), silent=T)
      if (exists ('a') ) {
        a<- as.data.frame (a)
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}

      if (exists ('a') == F ) {
        #in case there is no variation across the variable and then outliers function throws an error
        a <- data.frame (col1=rep (0,nrow(dat.environment)), col2 = rep (0,nrow(dat.environment) ))
        names (a) <- paste0(names(dat.environment)[i],c('-bxp','-rjk') )
        return (a)}



    })
    a.df <- do.call (cbind, a.df)
    choice.of.method.df <- a.df [,grep (names(a.df),pattern = outlier.method.biogeo)]

    if (class(.r.env)=="RasterLayer") {outlier.level <- choice.of.method.df}
    if (class(.r.env)=="RasterStack") {outlier.level <- round (rowMeans(choice.of.method.df) ,digits=2)}
    if (class(.r.env)=="RasterBrick") {outlier.level <- round (rowMeans(choice.of.method.df) ,digits=2)}

    # outlier.level.bxp <- round (rowMeans(a.df [,grep (names(a.df),pattern = '-bxp')]) ,digits=2)

    outlier.env <- (outlier.level >= .th.perc.outenv)*1
    comments.outlier.env <- paste (paste("Outlier",outlier.method.biogeo,".level="), outlier.level)

    out$envOutliers_bxp_test = outlier.env
    out$envOutliers_bxp_comments = comments.outlier.env

  }

  if (any (method %in% c('grubbs','all')))  {

    xydat <- df[,c(xf,yf)]
    spdf = sp::SpatialPoints(coords = xydat,proj4string =.projString )

    env.grubbs_tst = rep(0,times=nrow(out))
    env.grubbs.outl = occProfileR::findEnvOutliers(myPres = spdf,myEnv = .r.env,verbose = F)
    env.grubbs_tst [env.grubbs.outl] = 1
    env.grubbs_comment <- rep (paste('Pval= 1e-5'),times=nrow(out))

    out$envOutliers_Grubbs_test = env.grubbs_tst
    out$envOutliers_Grubbs_comments = env.grubbs_comment

  }


  out$envOutliers_score <- occProfileR:::.gimme.score (out)
  return (out)
}




#' @title Coordinate accuracy
#'
#' @description Detect environmental outliers using jacknife and boxplot
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param af the field where the geographic uncertainty is (in the same )
#' @param .r.env R enviroment
#' @param accept.threshold acceptance threshold for how much percentage of the Area of uncertainty in the cell we want to accept
#' @return
#' @keywords Analysis
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

geoEnvAccuracy <- function (df,
                         xf=x.field,yf=y.field,af=a.field,dsf=ds.field,ef =e.field,

                         method='all',

                         r.env,
                         accept.threshold.cell =0.5,
                         accept.threshold.env = 0.5,

                         bearing.classes=10,
                         distance.classes=5,
                         env.quantiles=c(0.3,0.7),

                         elev.threshold = 100,
                         raster.elevation=NULL,
                         verbose=F,
                         do=do.geoEnvAccuracy

                         # to be implemented : doParallel=T
                         ){



  #output results
  out = data.frame (geoenvLowAccuracy_lattice_test = NA,
                    geoenvLowAccuracy_lattice_comments = NA,
                    geoenvLowAccuracy_rasterCell_test = NA,
                    geoenvLowAccuracy_rasterCell_comments = NA,
                    geoenvLowAccuracy_envDiff_test =NA,
                    geoenvLowAccuracy_envDiff_comments =NA,
                    geoenvLowAccuracy_elevDiff_test =NA,
                    geoenvLowAccuracy_elevDiff_comments =NA,

                    geoenvLowAccuracy_score=NA
                    )[1:nrow (df),]

  row.names(x = out) <- NULL
  if (!do) {return (out)}

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


      cd_round_test = CoordinateCleaner::cd_round(x = xydat,lon = xf,lat = yf,ds=dsf,value='flagged',verbose=F,graphs=F)
      out$geoenvLowAccuracy_lattice_test = (!cd_round_test) * 1
      out$geoenvLowAccuracy_lattice_comments = ('Based on default values of CoordinateCleaner::cd_round')
  }

  #start methods requiring elevation
  #method elevation
  if (any(method %in% c('elevDiff','all'))) {

    #check and provide elevation
    if (is.null (ef) ) { out$elevDiff_test <- rep (NA,nrow(df))}

    if (!is.null(ef)){

      if (is.null (raster.elevation)) {
        if (verbose) print ('elevation raster not provided, downloading from SRTM')
        do.srtm.download = T
      } else {do.srtm.download=F}
      #if no elevation raster provided, download it
      if (do.srtm.download){
        raster.elevation = try  (expr = {.getSRTM (xydat=df,download=T, verbose=F)},silent = T)

      }

      #determining elevation test
      if (class(raster.elevation) == "RasterLayer" ) {
        elev.xy <- extract (raster.elevation, df[,c(xf,yf)])
        elev.difference <- abs (df[,ef] - elev.xy)
        out$elevDiff_test <- ifelse (elev.difference >= elev.threshold, 1, 0)
        out$elevDiff_comments <- paste ('Elev difference:', elev.difference)
      }
      #skipping elevation test
      if (class(raster.elevation) != "RasterLayer" ) {
        if(verbose) warning ("raster.elevation is not a raster; skipping elevDiff method")
        out$elevDiff_test <- (rep (NA,nrow(df)))
      }
    }

    }



  # THIS SHOULD BE THE LAST METHOD
  #check whether we need to do coordinate uncertainty methods because initial steps are time consuming
  if ( ! any (method %in% c('percDiffCell','envDeviation','all')) ) {
    #write final score
    out$geoenvLowAccuracy_score <-occProfileR:::.gimme.score (out)
    return (out)}
  if (is.null(af)) {
    print ('no coordinate accuracy/uncertainty field provided')
    #write final score
    out$geoenvLowAccuracy_score <-occProfileR:::.gimme.score (out)
    return (out)
    }

  #start methods requiring coordinate uncertainty
  if(length(af)>1) {df$new_accuracy = max(df[,af],na.rm=T); af = 'new_accuracy'}
  xydat = df[,c(xf,yf,af)]
  xydat$occCellids = cellFromXY(r.env,as.data.frame (df[,c(xf,yf)]))
  #get cell IDs of buffer points
  cellIds.directions = lapply (1:nrow (xydat), function (i){

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
      #perc.points.buffer = mean (points.buffer,na.rm=T)
      points.buffer.id
    })

    as.vector (out.bearing)
  })
  #start method percentage of Differnt cells around uncertainty
  if (any(method %in% c('percDiffCell','all'))){
    perc.diff <- sapply (1:nrow(xydat), function (i){
      pdiff <- mean ((cellIds.directions[[i]]==xydat$occCellids[i]) * 1, na.rm = T)
      pdiff
    })
    out$geoenvLowAccuracy_rasterCell_test = (perc.diff>accept.threshold.cell)*1
    out$geoenvLowAccuracy_rasterCell_comments = paste('Nbufferlocations=',bearing.classes*distance.classes)
  }
  #start method percentage of environmental differences around uncertainty
  if (any(method %in% c('envDeviation','all'))) {

    uniqueCellBuff <- unique (unlist (cellIds.directions))
    df.CellBuff = raster::extract(r.env,uniqueCellBuff,df=T)
    df.CellBuff$cellID = uniqueCellBuff

    targetCellUnique =   xydat$occCellids [!xydat$occCellids %in% uniqueCellBuff]
    df.CellTarget = raster::extract(r.env,targetCellUnique,df=T)
    df.CellTarget$cellID = targetCellUnique

    df.cells = rbind (df.CellBuff,df.CellTarget)

    env.test = lapply (1:nrow(xydat), function (i){


      NcellReps = table (cellIds.directions[[i]])

      df.env.cellBuffs  = lapply (1:length(NcellReps), function (m){
        cellID= as.numeric (names (NcellReps[m]))
        times = as.numeric (NcellReps[m])
        df.env.var = df.cells[which(df.cells$cellID ==cellID ), ]
        do.call("rbind", replicate(times, df.env.var, simplify = FALSE))

      })
      df.env.cellBuffs  = do.call(rbind,df.env.cellBuffs)
      df.env.cellBuffs  = df.env.cellBuffs[complete.cases(df.env.cellBuffs),]
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
      vars.outlier$all <- mean (as.numeric(vars.outlier[1,]),na.rm=T)
      vars.outlier

    })

    env.test = do.call (rbind, env.test)
    out$geoenvLowAccuracy_envDiff_test = (env.test$all > accept.threshold.env)*1
    out$geoenvLowAccuracy_envDiff_comments = paste('envQuantiles=',paste(env.quantiles,collapse = ','),';PercOutVariablesTh=',accept.threshold.env)

  }

  #write final score
  out$geoenvLowAccuracy_score <-occProfileR:::.gimme.score (out)
  return (out)
}

#'QUALITY GRADINGS
#'
#'Assess record quality
#' @param df Data.frame of species occurrences
#' @param qfield Field to add quality information in, default is "quality.comment"
#' @param new.comment New comment
#' @param separation.charachter Character to separate comments (?)
#' @return list
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
#'
#'
.add.to.qfield <- function (x,qfield='quality.comment',new.comment,separation.charachter=';'){
  stopifnot(is.data.frame(x))
  stopifnot(nrow(x) ==1)


  if(is.na (x[1,qfield])) {x[1,qfield] <- new.comment} else {x[1,qfield] <- paste (x[1,qfield],new.comment,collapse = separation.charachter)}

}



#' @title Centroid assessment using BIEN methods
#'
#' @descriptions Detect centroids in occurrences dataframe using BIEN methods
#' @details
#' @param 
#' @return list
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
# #' @export


centroid_assessment<-function(occurrences,centroid_data){
  
  occurrences$county<-as.character(occurrences$county)
  occurrences$state_province<-as.character(occurrences$state_province)
  occurrences$country<-as.character(occurrences$country)
  
  #message("This function assumes that occurrences are in wgs84 format, and then converts to whatever projection is used in the centroid data.  Only one projection at a time may be used in the centroid data, but multiple centroids per political division are fine.")
  
  #deal with nonsense points by dropping those lines
  
  if(length(  which(occurrences$latitude>90 | occurrences$latitude< -90) )>0  ){
    occurrences<-occurrences[-which(occurrences$latitude>90 | occurrences$latitude< -90),]}
  
  if(length(  which(occurrences$longitude>180 | occurrences$longitude< -180) )>0){
    occurrences<-occurrences[-which(occurrences$longitude>180 | occurrences$longitude< -180),]}
  
  if(any(is.na(occurrences$latitude) | (is.na(occurrences$longitude)))){
    occurrences<-occurrences[-which(is.na(occurrences$latitude) | (is.na(occurrences$longitude))),]  
    
  }
  
  temp_points<-SpatialPoints(coords = occurrences[c("longitude","latitude")],proj4string = CRS("+init=epsg:4326"))
  temp_points<-spTransform(x = temp_points,CRSobj = CRS(projargs = unique(as.character(centroid_data$projection))) )
  
  
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
      
      
      centroid_pt_i<-SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      
      
      country_data_i<-occurrences[which(occurrences$country==country_i),]
      
      new_data<-matrix(nrow = nrow(country_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      country_pts_i<-SpatialPointsDataFrame(coords = country_data_i[c("longitude","latitude")],data = as.data.frame(country_data_i$taxonobservation_id),proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- spDistsN1(pts = country_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(country_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-spDistsN1(pts = country_pts_i,pt = centroid_pt_i[p]) 
        
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
      
      centroid_pt_i<-SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      
      state_data_i<-occurrences[which(occurrences$country==state_i$country & occurrences$state_province==state_i$state_province),]
      
      new_data<-matrix(nrow = nrow(state_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      state_pts_i<-SpatialPointsDataFrame(coords = state_data_i[c("longitude","latitude")],data = as.data.frame(state_data_i$taxonobservation_id),proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- spDistsN1(pts = state_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(state_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-spDistsN1(pts = state_pts_i,pt = centroid_pt_i[p]) 
        
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
      
      
      centroid_pt_i<-SpatialPoints(coords = centroid_i[c("centroid_lon","centroid_lat")],proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      county_data_i<-occurrences[which(occurrences$country==county_i$country & occurrences$state_province==county_i$state_province & occurrences$county==county_i$county),]
      
      new_data<-matrix(nrow = nrow(county_data_i),ncol = nrow(centroid_i))
      new_data<-as.data.frame(new_data)
      
      county_pts_i<-SpatialPointsDataFrame(coords = county_data_i[c("longitude","latitude")],data = as.data.frame(county_data_i$taxonobservation_id),proj4string = CRS(projargs = unique(as.character(centroid_i$projection))))
      
      for(p in 1:nrow(centroid_i)){
        
        new_data[,p]<- spDistsN1(pts = county_pts_i,pt = centroid_pt_i[p]) / as.numeric(centroid_i$max_uncertainty[p] )
        
      }
      
      new_data_absolute<-matrix(nrow = nrow(county_data_i),ncol = nrow(centroid_i))
      new_data_absolute<-as.data.frame(new_data_absolute)
      for(p in 1:nrow(centroid_i)){
        
        new_data_absolute[,p]<-spDistsN1(pts = county_pts_i,pt = centroid_pt_i[p]) 
        
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





###### ========================================


#'Validates and clean occurrence data
#'
#'occurrence.profiler returns occurrence points withou the recommended metadata and attribution information.
#' @param sp.table data.frame. Object with the coordinate data.
#' @param sp.name character. Name of the species.
#' @param tableSettings list. Elements corresponding to different settings of the input occurrence table. 
#' @param analysisSettings list. Elements corresponding to different settings of the analysis functions . 
#' @param gradingSettings list. Elements corresponding to different settings of the analysis functions . 
#' @param writeoutSettings list. Elements corresponding to different settings of the analysis functions . 
#' @param r.env raster or rasterStack. Environmental data (e.g. typically climatic).
#' @param r.dem raster. Elevation data (in meters).
#' @return a list of two. First element is a dataframe with profiled occurrence records with their associated profiled labels. Second element is a dataframe with all outputs of the analysis implemented.
#' @note #There are several parameters in the function. The majority of them can be adjusted, but we also provide default values. We recommend those default values if the user is to use the geospatial data included in the package.
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export
occurrenceClassify <- function (
  sp.name,
  sp.table,
  r.env,
  
  tableSettings=NULL,
  analysisSettings=NULL,
  gradingSettings=NULL,
  writeoutSettings=NULL,
  r.dem=NULL,
  ntv.ctry=NULL,
  inv.ctry=NULL,
  resolveAlienCtry=F,
  resolveNativeCtry=F,
  interactiveMode=F,
  verbose = F){
  
  
  ########################################################################
  ### STEP 00: Load settings and study native and invasive countires
  ########################################################################
  
  defaults = system.file('ext/DefaultSettingsList.R'
                         ,package='occProfileR')
  source(defaults)
  
  #load table settings
  if (is.null(tableSettings)) { tableSettings = defaultSettings$tableSettings}
  taxonobservation.id = tableSettings$taxonobservation.id
  x.field             = tableSettings$x.field
  y.field             = tableSettings$y.field
  t.field             = tableSettings$t.field
  l.field             = tableSettings$l.field  
  c.field             = tableSettings$c.field
  e.field             = tableSettings$e.field
  a.field             =  tableSettings$a.field
  ds.field            =  tableSettings$ds.field
  
  #load analysisSettings
  if (is.null(analysisSettings)) { analysisSettings = defaultSettings$analysisSettings}
  
  coordinate.decimal.precision =  analysisSettings$geoSettings$coordinate.decimal.precision
  points.proj4string = analysisSettings$geoSettings$points.proj4string
  
  countries.shapefile = analysisSettings$rangeAnalysis$countries.shapefile
  countryfield.shapefile = analysisSettings$rangeAnalysis$countryfield.shapefile
  doRangeAnalysis = analysisSettings$rangeAnalysis$doRangeAnalysis
  excludeUnknownRanges = analysisSettings$rangeAnalysis$excludeUnknownRanges
  excludeNotmatchCountry = analysisSettings$rangeAnalysis$excludeNotmatchCountry
  doCountryRecordAnalysis = analysisSettings$rangeAnalysis$doCountryRecordAnalysis
  
  doCentroidDetection = analysisSettings$centroidAnalysis$doCentroidDetection
  methodCentroidDetection = analysisSettings$centroidAnalysis$methodCentroidDetection
  
  doHyperHumanDetection = analysisSettings$humanAnalysis$doHyperHumanDetection
  methodHyperHumanDetection = analysisSettings$humanAnalysis$methodHyperHumanDetection
  th.human.influence = analysisSettings$humanAnalysis$th.human.influence
  ras.hii = analysisSettings$humanAnalysis$ras.hii
  
  doInstitutionLocality = analysisSettings$institutionAnalysis$doInstitutionLocality
  methodInstitutionLocality = analysisSettings$institutionAnalysis$methodInstitutionLocality
  
  doGeoOutliers = analysisSettings$geooutliersAnalysis$doGeoOutliers
  methodGeoOutliers = analysisSettings$geooutliersAnalysis$methodGeoOutliers
  alpha.parameter = analysisSettings$geooutliersAnalysis$alpha.parameter
  
  doEnvOutliers = analysisSettings$envoutliersAnalysis$doEnvOutliers
  methodEnvOutliers = analysisSettings$envoutliersAnalysis$methodEnvOutliers
  th.perc.outenv= analysisSettings$envoutliersAnalysis$th.perc.outenv
  
  methodGeoEnvAccuracy=analysisSettings$accuracyAnalysis$methodGeoEnvAccuracy
  do.geoEnvAccuracy=analysisSettings$accuracyAnalysis$do.geoEnvAccuracy
  elev.quality.threshold=analysisSettings$accuracyAnalysis$elev.quality.threshold
  
  #load gradingSettings
  if (is.null(gradingSettings)) { gradingSettings = defaultSettings$gradingSettings}
  
  grading.test.type = gradingSettings$grading.test.type
  qualifiers = gradingSettings$qualifiers
  qualifier.label.scoping = gradingSettings$qualifier.label.scoping
  
  #load writeOUtSettings
  if (is.null(writeoutSettings)) { writeoutSettings = defaultSettings$writeoutSettings}
  output.dir = writeoutSettings$output.dir
  writeAllOutput = writeoutSettings$writeAllOutput
  write.simple.output =  writeoutSettings$write.simple.output
  write.full.output = writeoutSettings$write.full.output
  output.base.filename = writeoutSettings$output.base.filename
  
  #automatically resolve invasive and native countries for target species (not implemented yet) 
  
  if (interactiveMode & is.null (ntv.ctry) ){
    resolveNativeCtry <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species native range. Do you want to infer from global databases?")
  }
  if (interactiveMode & is.null (inv.ctry) ){
    resolveAlienCtry <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species alien range. Do you want to infer from global databases?")
  }
  
  if (any (resolveNativeCtry,resolveAlienCtry)){
    checkCountries <- occProfileR::nativeStatusCtry(spName = sp.name, resolveNative = resolveNativeCtry, resolveAlien = resolveAlienCtry ,verbose = verbose)
    if (resolveNativeCtry) {
      ntv.ctry <- c(ntv.ctry,checkCountries$ntvCtry)
      ntv.ctry <- unique (ntv.ctry)
    }
    if (resolveAlienCtry) {
      inv.ctry <- c(inv.ctry,checkCountries$invCtry)
      inv.ctry <- unique (inv.ctry)
    }
    
  }
  
  
  
  
  ########################################################################
  ### STEP 0: Initial checks
  ########################################################################
  
  
  if (missing (sp.table)) {stop('missing sp.table')}
  if (missing (r.env)) {stop('missing r.env')}
  
  
  
  ##############################################################################
  ### STEP 1: Data formatting and compatibility for biogeo and initial checks
  ##############################################################################
  
  #add fields necesary for initial table
  sp <- occProfileR:::.join.spname(sp.name)
  sp.table$Species <- sp
  
  sp.table2 <- occProfileR:::.checkfields(dat=sp.table,xf = x.field, yf=y.field,
                                          ef = e.field,tf = t.field,lf = l.field,
                                          cf = c.field,          idf = taxonobservation.id)
  
  dat <- occProfileR:::.addmainfields2(sp.table2,species = 'Species')
  dat$comments <- rep('', nrow(dat))
  
  #check data structure
  ck  <- occProfileR:::.checkdatastr2(dat,xf = x.field,yf=y.field)
  if (sum(ck$Present)!=10) {stop ("Error: required table fields could not be created")}
  
  
  
  #For development:
  #placeholder to convert other kind of coordinate data to numerical latlong
  #or whatever decided projection we want to work on
  
  #CHECK 4
  #ensure that all the geospatial data inputs are is in the same coordinate system
  potential.geosp.objects <- list(countries.shapefile,
                                  r.env,
                                  r.dem,
                                  ras.hii,
                                  points.proj4string)
  x <- sapply(potential.geosp.objects, is.null)
  actual.input.geosp.objects <- potential.geosp.objects[-x]
  actual.input.geosp.objects <- occProfileR:::.subsetlist.nonNULL(actual.input.geosp.objects)
  
  #the development should be in the direction of automatically check and transform
  #occProfileR:::.check.geospatial.data (list.geospatial.objects =actual.input.geosp.objects)
  
  #lapply (actual.input.geosp.objects, function (x) raster:::projection (x))
  
  
  #CHECK 5
  #get into precision of decimals
  for (i in c(x.field, y.field)){
    dat[,i] <- round (dat[,i], digits=coordinate.decimal.precision)
  }
  
  #CHECK 6
  #need tweaking: this will only happen if the data is in geographic coordinates, because the unites will be in minutes
  #res is is in minutes...therefore (30 arcsec = 0.5 min)
  #if this is a projected data, we need to redo stuff and also tweak funcitons in biogeo
  #maybe add a warning or something in check geospatial data? idk
  res.in.minutes <- res (r.env[[1]])[1] * 60
  
  ##############################################################################
  ### STEP 2: Quality H Filter : Identify records wo spatial info
  ##############################################################################
  Analysis.H <- filterMissing(df = dat,xf = x.field,yf = y.field)
  dat.Q.H <- Analysis.H$stay
  dat <- Analysis.H$continue
  
  #valid coordinates in geographic projections and zero zero issues anc decimal conversion issues
  if (as.character(points.proj4string) == 
      c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) {
    
    
    cc_val_test = CoordinateCleaner::cc_val(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=F)
    dat.Q.H1 = dat[!cc_val_test,]
    dat = dat[cc_val_test,]
    
    
    status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                              wfo = write.full.output,
                                                              wso = write.simple.output,
                                                              xf = x.field,
                                                              yf=y.field,
                                                              od = output.dir,
                                                              obf = output.base.filename,
                                                              sp=sp)
    
    #zero zero long lant
    cc_zero_test = CoordinateCleaner::cc_zero(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=F)
    dat.Q.H2 = dat[!cc_zero_test,]
    dat = dat[cc_zero_test,]
    
    status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                              wfo = write.full.output,
                                                              wso = write.simple.output,
                                                              xf = x.field,
                                                              yf=y.field,
                                                              od = output.dir,
                                                              obf = output.base.filename,
                                                              sp=sp)
    #error conversion of decimals 
    if(is.null(ds.field)) { 
      dat.tmp <- dat ;  dat.tmp$dataset = 'TemporaryDatasetName'
      cd_ddmm_test = CoordinateCleaner::cd_ddmm(x = dat.tmp,lon = x.field,lat = y.field,ds='dataset',value='flagged',verbose=F)
      dat.Q.H3 = dat[!cd_ddmm_test,]
      dat = dat[cd_ddmm_test,]
      rm (dat.tmp)
    }
    if(!is.null(ds.field)) { 
      cd_ddmm_test = CoordinateCleaner::cd_ddmm(x = dat.tmp,lon = x.field,lat = y.field,ds=ds.field,value='flagged',verbose=F)
      dat.Q.H3 = dat[!cd_ddmm_test,]
      dat = dat[cd_ddmm_test,]
      rm (dat.tmp)
    }
    
    status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                              wfo = write.full.output,
                                                              wso = write.simple.output,
                                                              xf = x.field,
                                                              yf=y.field,
                                                              od = output.dir,
                                                              obf = output.base.filename,
                                                              sp=sp)
  }
  
  
  #indicate issues of georeference and put them aside
  obj.issues = c('dat.Q.H','dat.Q.H1','dat.Q.H2','dat.Q.H3')
  obj.exist = sapply (obj.issues ,FUN = function (x) exists(x))
  obj.issues = obj.issues[obj.exist]
  if (any (obj.issues %in% c('dat.Q.H1','dat.Q.H2','dat.Q.H3'))) {
    if (length(obj.issues)>0){
      dat.excl.H  = lapply (obj.issues, function (x){
        
        
        if (nrow(get(x)) > 0 ) {
          a <- get (x)
          a$quality.grade <- 'H'
          a$Exclude <- 1 
          if (x == 'dat.Q.H')
            a$Reason <- "No coordinates"
          if (x == 'dat.Q.H1')
            a$Reason <- "No valid coords"
          if (x == 'dat.Q.H2')
            a$Reason <- "LonLat 0 0"
          if (x == 'dat.Q.H3')
            a$Reason <- "Coord decimal conversion error"
          return(a)
          
        } else { return (NULL)}
        
        
      })
      dat.excl.H = do.call(rbind, dat.excl.H)
    }
    if (exists('dat.excl.H') & (nrow (get('dat.excl.H')) > 0))
      dat.Q.H <- dat.excl.H
    if (exists ('dat.Q.H1')) {rm (dat.Q.H1)} 
    if (exists ('dat.Q.H2')) {rm (dat.Q.H2)} 
    if (exists ('dat.Q.H3')) {rm (dat.Q.H3)} 
    
    
  }
  
  #check outputs and escape ifneedbe //
  status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                            wfo = write.full.output,
                                                            wso = write.simple.output,
                                                            xf = x.field,
                                                            yf=y.field,
                                                            od = output.dir,
                                                            obf = output.base.filename,
                                                            sp=sp)
  
  if (is.list(status.out)) {return(status.out)}
  
  ##########################################################################################
  ### STEP 4: Filter Quality G : Identify duplicate records (in geographic space) to prevent pseudoreplicaton
  ##########################################################################################
  #indicate duplicates with exact coordinates
  Analysis.G <- duplicatesexcludeAnalysis(df = dat,
                                          xf = x.field,
                                          yf = y.field,
                                          resolution.in.minutes=res.in.minutes,
                                          raster.grid = r.env[[1]])
  
  #here we merge together the two kind of duplicates, but maybe we would like to keep the relative duplicates for other purposes
  dat.Q.G <- rbind (Analysis.G$Dups.Grid, Analysis.G$Dups.Exact)
  if (nrow(dat.Q.G) == 0 ) {rm (dat.Q.G)}
  if  (exists ('dat.Q.G')) { dat.Q.G$quality.grade <- 'G'}
  
  dat <- Analysis.G$continue
  
  #check outputs and escape ifneedbe //
  status.out=occProfileR:::.status.tracker.and.escaping(
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    obf = output.base.filename,
    sp=sp)
  
  
  if ( is.list  (status.out) ) {return (status.out)}
  
  
  ############################################################################
  ### STEP 5: SEA/TERRESTRIAL POTENTIAL REASSIGNMENT AND RECHECK DUPLICATES -- and potential for Quality G again(duplicates)
  ############################################################################
  #analysis of nearest cell next to the sea
  dat <- nearestcell3(dat=dat,rst = r.env, xf=x.field, yf=y.field)
  
  #check results and recheck dups ifneedbe
  if (class (dat) == 'list') {
    dat <- dat[[1]]
    moved.points <- dat[['moved']]
    if (!is.null(output.dir)){
      write.csv( moved.points,
                 paste0(output.dir,'/',sp,'_coastal_Reassignment.csv'))
    }
    
    
    
    ### RECHECK POTENTIAL DUPLICATES AGAIN AFTER REASSIGNATION
    Analysis.G.second.time=duplicatesexcludeAnalysis(df = dat,
                                                     xf = x.field,
                                                     yf = y.field,
                                                     resolution.in.minutes =
                                                       res.in.minutes)
    dat.Q.G.second.time <- rbind (Analysis.G.second.time$Dups.Grid,
                                  Analysis.G.second.time$Dups.Exact)
    if  (nrow (dat.Q.G.second.time) !=  0) {
      dat.Q.G.second.time$quality.grade <- 'G'}
    
    if (nrow (dat.Q.G.second.time) !=  0) {
      dat.Q.G <- rbind(dat.Q.G,dat.Q.G.second.time)}
    
    rm (dat.Q.G.second.time)
  }
  
  #check outputs and escape if need be //
  status.out <- occProfileR:::.status.tracker.and.escaping(
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    obf = output.base.filename,
    sp=sp)
  
  
  if ( is.list(status.out) ) {return(status.out)}
  
  
  ##########################################################################
  ### STEP 6: Filter Quality F Country selection
  
  Analysis.F <-occProfileR:::countryStatusRangeAnalysis(df=dat,
                                                        xf = x.field,
                                                        yf = y.field,
                                                        .ntv.ctry = ntv.ctry,
                                                        .inv.ctry = inv.ctry,
                                                        .c.field = c.field,
                                                        .countries.shapefile =
                                                          countries.shapefile,
                                                        cfsf = countryfield.shapefile,
                                                        .points.proj4string =
                                                          points.proj4string,
                                                        excludeUnknownRanges = excludeUnknownRanges,
                                                        excludeNotmatchCountry = excludeNotmatchCountry,
                                                        doRangeAnalysis = doRangeAnalysis,
                                                        verbose = F
  )
  
  if  (nrow (Analysis.F$stay) != 0 ) {   dat.Q.F <- Analysis.F$stay ; dat.Q.F$quality.grade <- 'F'}
  dat <- Analysis.F$continue
  
  #check outputs and escape ifneedbe //
  status.out <- occProfileR:::.status.tracker.and.escaping (
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    obf = output.base.filename,
    sp=sp)
  
  
  if ( is.list  (status.out) ) {return (status.out)}
  
  ############################################################################
  ### STEP 7:  Quality A-E Environmental and Geographical outliers  - analysis chunk
  ############################################################################
  
  #ANALYSIS ELEMENTS
  #this is important for development, need to specify the number of ELEMENTS of analysis
  #to sumarize the results later on we will need that number
  N_analysis = 8
  
  
  
  
  ### ELEMENT 1: CENTROID ISSUE DETECTION
  Analysis.1     <- centroidDetection (.r.env = r.env,
                                       df = dat,
                                       xf = x.field,
                                       yf = y.field,
                                       cf = c.field,
                                       idf = taxonobservation.id,
                                       .ntv.ctry = ntv.ctry,
                                       .inv.ctry = inv.ctry,
                                       .points.proj4string =
                                         points.proj4string,
                                       .countries.shapefile = countries.shapefile,
                                       method = methodCentroidDetection,
                                       do= doCentroidDetection)
  
  ### ELEMENT 2: HYPER-HUMAN ENVIRONMENT
  Analysis.2 <- hyperHumanDetection(df = dat,
                                    xf = x.field,
                                    yf = y.field,
                                    .points.proj4string =points.proj4string,
                                    ras.hii = ras.hii,
                                    .th.human.influence =th.human.influence,
                                    do = doHyperHumanDetection)
  
  
  
  
  
  ### ELEMENT 3: BOTANICAL GARDEN PLACEMENT -- FROM LOCALITY NAME
  Analysis.3 <- institutionLocality (df=dat,lf=l.field,xf = x.field,yf=y.field,
                                     do = doInstitutionLocality,
                                     method = methodInstitutionLocality)
  
  ### ELEMENT 4: GEOGRAPHICAL OUTLIER
  Analysis.4 <- geoOutliers(df=dat,
                            xf=x.field,
                            yf=y.field,
                            .alpha.parameter = alpha.parameter,
                            #do=doGeoOutliers,
                            method = methodGeoOutliers,
                            .projString = points.proj4string)
  
  
  ### ELEMENT 5: ADD the geoIndicator of wrong country if that has not been determined to be an exclusive dimension
  
  if (!excludeUnknownRanges & doRangeAnalysis){
    
    a = dat$countryStatusRangeAnalysis_wrongNTV_test ; a[is.na(a)] <- 0
    b = dat$countryStatusRangeAnalysis_wrongINV_test ; b[is.na(b)] <- 0
    
    Analysis.5 <-  a+b
    Analysis.5 <-  ifelse(test = Analysis.5>1,yes = 1,no = 0)
    Analysis.5 <- data.frame (countryStatusRangeAnalysis_notinKnownRange_test=Analysis.5, countryStatusRangeAnalysis_notinKnownRange_comments='')
    Analysis.5$unknownRange_score <- occProfileR:::.gimme.score(Analysis.5)
    row.names(Analysis.5) <- NULL
    rm(a);rm(b)
  }
  if (excludeUnknownRanges | !doRangeAnalysis){
    Analysis.5 <- data.frame (countryStatusRangeAnalysis_notinKnownRange=NA,
                              countryStatusRangeAnalysis_notinKnownRange_comments='unkonown ranges are classifed as F directly',
                              unknownRange_score=NA)[1:nrow(dat),]
    row.names(Analysis.5) <- NULL
  }
  
  ### ELEMENT 6: ADD the geoIndicator of country recorded different than country coordinates
  if (!excludeNotmatchCountry & doCountryRecordAnalysis){
    Analysis.6 = dat$countryStatusRangeAnalysis_wrongCTRY_test
    Analysis.6 <- data.frame (countryStatusRangeAnalysis_wrongReportCtry_test=Analysis.6,
                              countryStatusRangeAnalysis_wrongReportCtry_comments='')
    Analysis.6$wrongReportCtry_score <- occProfileR:::.gimme.score(Analysis.6)
  }
  if (excludeNotmatchCountry | !doCountryRecordAnalysis){
    Analysis.6 <- data.frame (countryStatusRangeAnalysis_wrongReportCtry=NA,
                              countryStatusRangeAnalysis_wrongReportCtry_comments=NA,
                              wrongReportCtry_score= NA)[1:nrow(dat),]
    row.names(Analysis.6) <- NULL
  }
  
  ### ELEMENT 7: ENVIRONMENTAL OUTLIER
  
  Analysis.7 <- envOutliers (.r.env=r.env,
                             df= dat, xf=x.field,
                             yf =y.field,
                             .th.perc.outenv = th.perc.outenv,
                             .sp.name = sp.name,
                             .projString = points.proj4string,
                             method = methodEnvOutliers,
                             do = doEnvOutliers)
  
  
  ### ELEMENT 8: Coordinate accuracy
  Analysis.8 <- geoEnvAccuracy(df=dat,
                               xf = x.field,
                               yf = y.field,
                               af = a.field,
                               dsf= ds.field,
                               r.env = r.env,
                               ef= e.field,
                               raster.elevation = r.dem,
                               do = do.geoEnvAccuracy,method = methodGeoEnvAccuracy)
  
  
  
  
  
  
  
  
  ### SUMMARY ANALYSIS RESULTS
  list.analysis = list()
  for (i in 1:N_analysis){
    list.analysis[[i]] <- get (paste0('Analysis.',i))
  }
  
  df.qualityAssessment <- do.call(cbind,list.analysis)
  row.names(df.qualityAssessment) <- NULL
  
  
  
  ############################################################################
  #### STEP 8:  Quality A-E Environmental and Geographical outliers
  # assignation of grade to a record
  ###########################################################################
  if (gradingSettings$grading.test.type == 'strict')   {test.strictness.value = 0 }
  if (gradingSettings$grading.test.type == 'majority') {test.strictness.value = 0.5 }
  if (gradingSettings$grading.test.type == 'relaxed')  {test.strictness.value = 0.6 }
  #start assigning values depending on conditions
  suppressMessages (attach(df.qualityAssessment))
  df.qualityAssessment$quality.grade <- NA
  
  #load lazylogic function
  .lazylogic <- function (e,...){
    
    o <- eval(parse(text=e))
    ifelse(is.na(o),F,o)
    
  }
  
  grade.E <- (.lazylogic (e = 'hyperHumanDetection_score >= test.strictness.value')) | (.lazylogic (e = 'institutionLocality_score >= test.strictness.value')) | (.lazylogic (e = 'centroidDetection_score >= test.strictness.value')) | (.lazylogic (e = 'unknownRange_score >= test.strictness.value')) | (.lazylogic (e = 'wrongReportCtry_score >= test.strictness.value') | (.lazylogic (e = 'envOutliers_missingEnv_score >= test.strictness.value')))
  df.qualityAssessment$quality.grade[grade.E] <- 'E'
  to.continue <- is.na(df.qualityAssessment$quality.grade)
  
  grade.D <- (.lazylogic (e = 'geoOutliers_score >= test.strictness.value') ) &  (.lazylogic ('envOutliers_score  >= test.strictness.value') )
  df.qualityAssessment$quality.grade[to.continue & grade.D] <- 'D'
  to.continue <- is.na(df.qualityAssessment$quality.grade)
  
  grade.C <- .lazylogic ('envOutliers_score  >= test.strictness.value')
  df.qualityAssessment$quality.grade[to.continue & grade.C] <- 'C'
  to.continue <- is.na(df.qualityAssessment$quality.grade)
  
  grade.B <- .lazylogic (e = 'geoOutliers_score >= test.strictness.value') | (.lazylogic (e = 'geoenvLowAccuracy_score >= test.strictness.value') ) 
  df.qualityAssessment$quality.grade[to.continue & grade.B] <- 'B'
  to.continue <- is.na(df.qualityAssessment$quality.grade)
  
  df.qualityAssessment$quality.grade[to.continue] <- 'A'
  
  detach(df.qualityAssessment)
  
  #add labels to the data
  dat <- cbind(dat,df.qualityAssessment)
  
  
  
  ############################################################################
  ### STEP 9: BUILD FULL dataframe
  ############################################################################
  previous.Q.objects <- grep(pattern = 'dat.Q.',ls(),value = T)
  full.qaqc <- dat
  for (o in previous.Q.objects){
    full.qaqc <- plyr::rbind.fill(full.qaqc, get(o))
  }
  
  
  ###########################################################################
  ### STEP 10: GET QUALIFIERS (attributes you want to add to your label that help decide its use)
  ###########################################################################
  #check if user does not want the qualifiers, and return output
  if (gradingSettings$qualifiers == F) {
    
    full.qaqc$qualifiers <- NA
    full.qaqc$quality.label <- full.qaqc$quality.grade
    
    #WRITE OUTPUTS AND END THE PROCESS
    if (write.full.output==T) {
      write.csv (full.qaqc,   paste0(output.dir,'/',sp,'_long_' ,
                                     output.base.filename,'.csv'),row.names = F)
    }
    
    short.qaqc<-full.qaqc[,c(x.field,y.field,'quality.grade',
                             'qualifiers','quality.label')]
    if (write.simple.output==T) {
      write.csv (short.qaqc,
                 paste0(output.dir,'/',sp,'_short_' ,output.base.filename,
                        '.csv'),row.names = F)
    }
    
    output.function <- list (occ_full_profile=full.qaqc,
                             occ_short_profile=short.qaqc)
    return (output.function)
  }
  
  
  
  
  
  #get the qualifiers
  tag1 <- occProfileR:::.qualifier.invasive.range(df=full.qaqc,
                                                  tag = 'i',
                                                  qualifier.lab.scp = qualifier.label.scoping)
  tag2 <- occProfileR:::.qualifier.native.range(df=full.qaqc,
                                                tag = 'n',
                                                qualifier.lab.scp = qualifier.label.scoping)
  
  #we will also add if we have a time stamp
  tag3 <-  occProfileR:::.qualifier.timestamp(df=full.qaqc,
                                              tf = t.field,
                                              tag = 'T',
                                              qualifier.lab.scp = qualifier.label.scoping)
  
  #we will also add if we have the right precision
  tag4 <- occProfileR:::.qualifier.precision(df=full.qaqc,tag='p',
                                             .s = res.in.minutes,
                                             .e = res.in.minutes,
                                             xf = x.field,
                                             yf = y.field,
                                             qualifier.lab.scp = qualifier.label.scoping)
  
  #we will also add an elevation threshold level
  tag5 <- occProfileR:::.qualifier.elevation (df=full.qaqc,tag='e',
                                              qualifier.lab.scp = qualifier.label.scoping 
  )
  
  #get all tags together
  tags.applied <- grep('tag',ls(),value=T)
  tag.output <- lapply (1:length(tag1), function (i){
    
    row.tag <- unlist(lapply (tags.applied, function (x){a <- get(x)[i] }))
    row.tag <- paste(row.tag,collapse = '-')
    
    return(row.tag)
    
  })
  tag.output <- unlist (tag.output)
  tag.output <- lapply (tag.output, function (x){
    
    row.tag <- strsplit (x,split='-')[[1]]
    row.tag <- row.tag [row.tag!="NA"]
    if (length(row.tag)==0){return(NA)}
    if (length(row.tag)==1){return(row.tag)}
    if (length(row.tag)>1) {row.tag <- paste(row.tag,collapse = '-'); return (row.tag)}
    
  })
  tag.output <- unlist (tag.output)
  full.qaqc$qualifiers <- tag.output
  
  #collapse to a unique label
  full.qaqc$quality.label <- occProfileR:::.paste3 (as.character(full.qaqc$quality.grade),full.qaqc$qualifiers,sep = '/')
  
  
  ########################################################################
  ### STEP 11: WRITE THE OUTPUTS
  ########################################################################
  if (write.full.output==T) {
    write.csv (full.qaqc,  paste0(output.dir,'/',
                                  output.base.filename,'_',sp,'_long.csv'),
               row.names = F)
  }
  
  short.qaqc<-full.qaqc[,c(x.field,y.field,'quality.grade',
                           'qualifiers','quality.label',
                           "countryStatusRangeAnalysis_score", "centroidDetection_score" ,"hyperHumanDetection_score" ,"institutionLocality_score",
                           "geoOutliers_score","unknownRange_score", "wrongReportCtry_score","envOutliers_score" ,"geoenvLowAccuracy_score")]
  if (write.simple.output==T) {
    write.csv (short.qaqc,
               paste0(output.dir,'/',output.base.filename,'_',sp,
                      '_short.csv'),row.names = F)
  }
  
  output.function <- list (occ_full_profile=full.qaqc, occ_short_profile=short.qaqc)
  return (output.function)
  
  
}

