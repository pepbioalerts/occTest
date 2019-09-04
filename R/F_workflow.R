
###### ========================================


#'Validates and clean occurrence data
#'
#'occurrence.profiler returns occurrence points withou the recommended metadata and attribution information.
#' @param output.dir Character. The folder where the output will be dumped.
#' @param sp.table data.frame. Object with the coordinate data.
#' @param sp.name character. Name of the species.
#' @param taxonobservation.id characther. field with the taxon observation identifier (default = NULL).
#' @param r.env raster or rasterStack. Environmental data (e.g. typically climatic).
#' @param r.dem raster. Elevation data (in meters).
#' @param countries.shapefile SpatialPolygonsDataframe. Country borders with associated data on country name, etc. etc.
#' @param countryfield.shapefile character. Name of the field in **_countries.shapefile_** where the 3 letter ISO is stored. Default is "ISO".
#' @param ntv.ctry character. Three letter code of the country/countries were the species is native. If set to NULL then all occurrences will be considered as in the native range.
#' @param inv.ctry character. Three letter code of the country/countries were the species is invasive If set to NULL then the species will never considered invasive anywhere.
#' @param ras.hii raster. Raster with associated values informing on the degree of anthropogenic influence in a place.
#' @param x.field character. Name of the field where the x coordinate is stored (typically longitude).
#' @param y.field character. Name of the field where the y coordinate is stored (typically latitude).
#' @param t.field character. Name of the field where the date of data collection is stored in the original dataset. Default NULL.
#' @param c.field character. Name of the field where the registred country of data collection is stored in the original dataset. Default NULL.
#' @param l.field character. Name of the field where the toponim/location of data collection is stored in the original dataset. Default NULL.
#' @param e.field character. Name of the field where the elevation of data collection is stored in the original dataset. Default NULL.
#' @param coordinate.decimal.precision numeric. Decimal precision of the coordinates. Default is 4.
#' @param points.proj4string CRS object. Projection of the orginal presence coordinates. Set toCRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"). This option can't be changed at the moment.
#' @param alpha.parameter numeric. Parameter for determining the alpha.hull for geographic outlier detection. Default is 2.
#' @param th.human.influence numeric. Maximum level of anthropogenic influence for anthropogenic issues detection. Default is 45.
#' @param th.perc.outenv numeric. Ratio of variables considered outliers for evnironmental spaces issues detection. Default is 0.2.
#' @param elev.quality.threshold numeric. Maximum difference between recorded and expected elevation allowed for elevation quality assessement. Default is 100.
#' @param qualifiers logical. Compute qualifiers. Deafault is T.
#' @param qualifier.label.scoping character. Quality grades over which to compute the qualifiers. Default is "A".
#' @param write.simple.output logical. Write output table with basic information of the grading. Deafault is F.
#' @param write.full.output logical. Write output table with full information of the grading. Deafault is F.
#' @param output.base.filename character. Base filenames for outputs. Meaningfull if write outputs are set to T. Default is "QAQC"
#' @return a list of two. First element is a dataframe with profiled occurrence records with their associated profiled labels. Second element is a dataframe with all outputs of the analysis implemented.
#' @note #There are several parameters in the function. The majority of them can be adjusted, but we also provide default values. We recommend those default values if the user is to use the geospatial data included in the package.
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export
occurrenceProfile <- function (
          output.dir=NULL,
          
          sp.table,
          sp.name,
          
          taxonobservation.id = NULL,
          x.field = 'x',
          y.field = 'y',
          t.field = NULL,
          l.field = NULL,
          c.field = NULL,
          e.field = NULL,
          a.field = NULL, #coordinate uncertainty in meters
          ds.field = NULL, #dataset field identifier

          r.env,
          r.dem=NULL,
          countries.shapefile=NULL,
          countryfield.shapefile = 'ISO',
          
          ntv.ctry=NULL,
          inv.ctry=NULL,
          ras.hii=NULL,
          coordinate.decimal.precision = 4,
          points.proj4string=sp:::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),

          #analysis to perform
          doRangeAnalysis=T,
          excludeUnknownRanges= F,
          excludeNotmatchCountry= F,
          doCountryRecordAnalysis=T,
          doCentroidDetection=T,
          methodCentroidDetection='all',
          doHyperHumanDetection=T,
          methodHyperHumanDetection='all',
          doInstitutionLocality=T,
          methodInstitutionLocality='all',
          doGeoOutliers=T,
          methodGeoOutliers='all',
          doEnvOutliers=T,
          methodEnvOutliers='all',
          methodGeoEnvAccuracy='all',
          do.geoEnvAccuracy=T,

          #detail.parameters
          alpha.parameter = 2,
          th.human.influence = 45,
          th.perc.outenv =  0.2,
          elev.quality.threshold = 100,

          #quality assessment parameters
          grading.test.type = 'majority', #other options are 'strict' 'relaxed'
          #test.strictness.value = 0.5,
          qualifiers=T,
          qualifier.label.scoping=c('A','B','C','D','E'),

          #writing outputs
          write=F, #overwrites write.simple.output, write.full.output
          write.simple.output=F,
          write.full.output=F,
          output.base.filename="QAQC",
          verbose = F){

  #functions params to rerun
  # output.dir; sp.table; sp.name;r.env; r.dem;
  # countries.shapefile;countryfield.shapefile = 'ISO';ntv.ctry; inv.ctry; ras.hii;
  # taxonobservation.id = 'ID_GRAFIC';x.field = 'x'; y.field = 'y'; t.field = NULL;l.field = NULL;c.field =NULL;e.field = NULL;
  # coordinate.decimal.precision = 4;points.proj4string=sp:::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
  # alpha.parameter = 2;th.human.influence = 45;th.perc.outenv =  0.2;elev.quality.threshold = 100;
  # write.simple.output=T
  # write.full.output=T
  # excludeUnknownRanges= F ; excludeNotmatchCountry= F ; doCentroidDetection=T ; methodCentroidDetection='all'
  # doHyperHumanDetection=T ; doInstitutionLocality=T ; methodInstitutionLocality='all'; doGeoOutliers=T ; methodGeoOutliers='all'
  # doRangeAnalysis = T; doCountryRecordAnalysis = T ; doEnvOutliers=T
  # output.base.filename="QAQC"
  # grading.test.type='majority'
  # methodEnvOutliers ='all'



  ########################################################################
  ### STEP 0: Initial checks
  ########################################################################
  #stopifnot(is.character (output.base.filename))
  #if(dir.exists(paths = output.dir) == F) {dir.create (output.dir,
  #                                                     showWarnings = T,
  #                                                     recursive = T)}
  #stopifnot(dir.exists(paths = output.dir) == T)



  ##############################################################################
  ### STEP 1: Data formatting and compatibility for biogeo and initial checks
  ##############################################################################

  #CHECK 0
  #all necessary data is there

  #if (missing (output.dir)) {stop('missing output.dir')}

  #if (doRangeAnalysis) {
    #if (is.null (countries.shapefile)) {stop('missing countries.shapefile')}
    #if (is.null (ntv.ctry)) {stop('missing ntv.ctry')}
    #if (is.null (inv.ctry)) {stop('missing inv.ctry')}
  #}
  #if () {if (missing (r.dem)) {stop('missing r.dem')}}
  

  if (missing (sp.table)) {stop('missing sp.table')}
  if (missing (sp.name)) {stop('missing sp.name')}
  if (missing (r.env)) {stop('missing r.env')}
  if (doHyperHumanDetection) {if (is.null (ras.hii)) {stop('missing ras.hii')}}
  if (is.null(countries.shapefile)){
    #try to have the countries downloaded the countries downloaded from an external point
    
    }

  if (!write) {write.full.output<-F;write.simple.output<-F}

  #CHECK 1
  #add fields necesary for initial table
  sp <- occProfileR:::.join.spname(sp.name)
  sp.table$Species <- sp

  sp.table2 <- occProfileR:::.checkfields(dat=sp.table,xf = x.field, yf=y.field,
                            ef = e.field,tf = t.field,lf = l.field,
                            cf = c.field,          idf = taxonobservation.id)
  #if (is.null (e.field) )         {e.field<- 'elev' }
  #if (is.null (l.field) )         {l.field<- 'locality'}
  #if (is.null (c.field) )         {c.field<- 'countryRecorded' }
  #if (is.null (t.field) )         {t.field<- 'time'}

  dat <- occProfileR:::.addmainfields2(sp.table2,species = 'Species')
  dat$comments <- rep('', nrow(dat))

  #CHECK 2
  #check data structure
  ck  <- occProfileR:::.checkdatastr2(dat,xf = x.field,yf=y.field)
  if (sum(ck$Present)!=10) {stop ("Error: required table fields could not be created")}


  #CHECK 3
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
                                       .ntv.ctry = ntv.ctry,
                                       .inv.ctry = inv.ctry,
                                       .points.proj4string =
                                         points.proj4string,
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
                            do=doGeoOutliers,
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
  if (grading.test.type == 'strict')   {test.strictness.value = 0 }
  if (grading.test.type == 'majority') {test.strictness.value = 0.5 }
  if (grading.test.type == 'relaxed')  {test.strictness.value = 0.6 }
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
  if (qualifiers == F) {

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
  defaultSettings = occProfileR::defaultSettings()
  
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


