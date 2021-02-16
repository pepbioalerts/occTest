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
occurrenceTests <- function (
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
  outPath=NULL,
  verbose = F,
  doParallel=F,
  mc.cores=2){
  
  #set timer
  tictoc:::tic()
  
  ########################################################################
  ### STEP 00: Initial checks
  ########################################################################
  #set timer
  tictoc:::tic('Initial checks and formatting')
  print ('Initial checks and formatting started...')
  
  if (missing (sp.table)) {stop('missing sp.table')}
  if (missing (r.env)) {stop('missing r.env')}
  if (! pingr:::is_online()) { stop('You seem not to have Internet connection. This package requires internet connection for several tests. Please go online')}
  
  
  ########################################################################
  ### STEP 0: Load settings and study native and invasive countries
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
  
  doHumanDetection = analysisSettings$humanAnalysis$doHumanDetection
  methodHumanDetection = analysisSettings$humanAnalysis$methodHumanDetection
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
  
  #load writeOutSettings
  if (is.null(writeoutSettings)) { writeoutSettings = defaultSettings$writeoutSettings}
  output.dir = writeoutSettings$output.dir
  if (is.null(output.dir)) {output.dir = getwd() ; message('Output directory set to your current working directory')}
  writeAllOutput = writeoutSettings$writeAllOutput
  write.simple.output =  writeoutSettings$write.simple.output
  write.full.output = writeoutSettings$write.full.output
  output.base.filename = writeoutSettings$output.base.filename
  
  
  ##############################################################################
  ### STEP 1a: Data formatting and compatibility for biogeo and initial checks
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
  
  
  tictoc:::toc()
  
  ##############################################################################
  #### NOT IMPLEMENTED YET !!!!!!!!! TO DEVELOP WITH BIEN GNRS/NRS ,...
  ### STEP 1b [OPTIONAL]: Automatically solve native or invasive range 
  ##############################################################################
  
  #whatever the want we are going to set it to False
  interactiveMode= F
  resolveNativeCtry = F
  resolveAlienCtry = F
  #set timer
  #tictoc:::tic('Resolve native and invasive countries')
  
  #automatically resolve invasive and native countries for target species (not implemented yet) 
  
  if (interactiveMode & is.null (ntv.ctry) ){
    resolveNativeCtry <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species native range. Do you want to infer from global databases?")
  }
  if (interactiveMode & is.null (inv.ctry) ){
    resolveAlienCtry <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species alien range. Do you want to infer from global databases?")
  }
  if (any (resolveNativeCtry,resolveAlienCtry)){
    if (verbose) {print ("**** RESOLIVNG NATIVE AND INVASIVE RANGES ****")}
    
    xydatTemporary = dat[,c(x.field,y.field)]
    xydatTemporary = xydatTemporary[complete.cases(xydatTemporary),]
    
    
    checkCountries <- occProfileR:::nativeStatusCtry(spName = sp.name, xydat=xydatTemporary, resolveNative = resolveNativeCtry, resolveAlien = resolveAlienCtry ,verbose = verbose)
    if (resolveNativeCtry) {
      ntv.ctry <- c(ntv.ctry,checkCountries$ntvCtry)
      ntv.ctry <- unique (ntv.ctry)
    }
    if (resolveAlienCtry) {
      inv.ctry <- c(inv.ctry,checkCountries$invCtry)
      inv.ctry <- unique (inv.ctry)
    }
    
  }
  #tictoc:::toc()
  
  ##############################################################################
  ### STEP 2: Quality H Filter : Identify records wo spatial info
  ##############################################################################
  #set timer
  tictoc:::tic('Filter major coordinate Issues')
  print('Filter major coordinate Issues statrted...')
  
  if (verbose) {print ("**** RESOLVING QUALITY FILTER: records wo Spatial Info ****")}
  
  Analysis.H <- filterMissing(df = dat,xf = x.field,yf = y.field)
  dat.Q.H <- Analysis.H$stay
  dat <- Analysis.H$continue
  
  #valid coordinates in geographic projections and zero zero issues anc decimal conversion issues
  if (as.character(points.proj4string) %in% 
      c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',"+proj=longlat +datum=WGS84 +no_defs") ) {
    
    
    coordIssues_invalidCoord_value =! CoordinateCleaner::cc_val(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=F)
    dat$coordIssues_invalidCoord_value = coordIssues_invalidCoord_value
    coordIssues_invalidCoord_test = as.logical (coordIssues_invalidCoord_value)
    dat$coordIssues_invalidCoord_test = coordIssues_invalidCoord_test
    
    dat.Q.H1 = dat[coordIssues_invalidCoord_test,]
    dat = dat[!coordIssues_invalidCoord_test,]
    
    
    status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                              wfo = write.full.output,
                                                              wso = write.simple.output,
                                                              xf = x.field,
                                                              yf=y.field,
                                                              od = output.dir,
                                                              obf = output.base.filename,
                                                              sp=sp)
    
    #zero zero long lant
    coordIssues_zeroCoord_value = ! CoordinateCleaner::cc_zero(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=F)
    dat$coordIssues_zeroCoord_value = coordIssues_zeroCoord_value
    coordIssues_zeroCoord_test = as.logical (coordIssues_zeroCoord_value)
    dat$coordIssues_zeroCoord_test = coordIssues_zeroCoord_test
    dat.Q.H2 = dat[coordIssues_zeroCoord_test,]
    dat = dat[!coordIssues_zeroCoord_test,]
    
    status.out <- occProfileR:::.status.tracker.and.escaping (dataset.to.continue = dat,
                                                              wfo = write.full.output,
                                                              wso = write.simple.output,
                                                              xf = x.field,
                                                              yf=y.field,
                                                              od = output.dir,
                                                              obf = output.base.filename,
                                                              sp=sp)
    #error conversion of decimals 
    dat.tmp <- dat 
    if(is.null(ds.field)) { warning('No dataset field provided, considering everything as a unique dataset')
      dat.tmp$MyInventedCommonDataset = 'TemporaryDatasetName' 
      ds.field = 'MyInventedCommonDataset'}
    coordIssues_coordConv_value = ! CoordinateCleaner::cd_ddmm(x = dat.tmp,lon = x.field,lat = y.field,ds=ds.field , value='flagged',verbose=F)
    dat$coordIssues_coordConv_value = coordIssues_coordConv_value
    coordIssues_coordConv_test = as.logical (coordIssues_coordConv_value)
    dat$coordIssues_coordConv_test = coordIssues_coordConv_test
    
    dat.Q.H3 = dat[coordIssues_coordConv_test,]
    dat = dat[!coordIssues_coordConv_test,]
    if (ds.field == 'MyInventedCommonDataset') ds.field <- NULL
    rm (dat.tmp)
    
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
          #a$quality.grade <- 'H'
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
    if (exists('dat.excl.H') & (!is.null (get('dat.excl.H')))) {
      if (nrow (dat.excl.H)>0) dat.Q.H <- dat.excl.H
      rm (dat.excl.H)
    }
    if (exists('dat.excl.H') & (is.null (get('dat.excl.H')))) {
      rm (dat.excl.H)
    }
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
  tictoc:::toc()
  ##########################################################################################
  ### STEP 4: Filter Quality G : Identify duplicate records (in geographic space) to prevent pseudoreplicaton
  ##########################################################################################
  #set timer
  tictoc:::tic('Filter duplicates')
  print('Filter duplicates started')
  if (verbose) {print ("**** RESOLVING: duplicates ****")}
  
  #indicate duplicates with exact coordinates
  Analysis.G <- duplicatesexcludeAnalysis(df = dat,
                                          xf = x.field,
                                          yf = y.field,
                                          resolution.in.minutes=res.in.minutes,
                                          raster.grid = r.env[[1]])
  
  #here we merge together the two kind of duplicates, but maybe we would like to keep the relative duplicates for other purposes
  dat.Q.G <- dplyr::bind_rows (Analysis.G$Dups.Grid, Analysis.G$Dups.Exact)
  if (nrow(dat.Q.G) == 0 ) {rm (dat.Q.G)}
  #if  (exists ('dat.Q.G')) { dat.Q.G$quality.grade <- 'G'}
  
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
  
  tictoc:::toc()
  ############################################################################
  ### STEP 5: SEA/TERRESTRIAL POTENTIAL REASSIGNMENT AND RECHECK DUPLICATES -- and potential for Quality G again(duplicates)
  ############################################################################
  #set timer
  tictoc:::tic('Resolving coastal reassignment')
  print ('Resolving coastal reassignment started...')
  if (verbose) {print ("**** RESOLVING : sea/terrestrial reassignment ****")}
  #analysis of nearest cell next to the sea
  dat <- nearestcell3(dat=dat,rst = r.env, xf=x.field, yf=y.field)
  
  #check results and recheck dups ifneedbe
  if (class (dat) == 'list') {
    dat <- dat[[1]]
    moved.points <- dat[['moved']]
    if (!is.null(output.dir)){
      sp.name2= occProfileR:::.join.spname(sp.name)
      odir = paste0(output.dir,'/',sp.name2)
      dir.create(odir,showWarnings = F,recursive = T)
      write.csv( moved.points,
                 paste0(odir,'/',sp.name2,'_coastal_Reassignment.csv'))
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
      #dat.Q.G.second.time$quality.grade = 'G'
      dat.Q.G = rbind(dat.Q.G,dat.Q.G.second.time)
      dat = Analysis.G.second.time$continue
    }
    
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
  tictoc:::toc()
  
  
  ##########################################################################
  ### STEP 6: Filter Quality F Country selection
  
  #set timer
  tictoc:::tic('Resolving countryStatusRange Analysis')
  print('Resolving countryStatusRange Analysis started...')
  if (verbose) {print ("**** RESOLVING QUALITY FILTER F: CountryStatus ranege analysis (invasive?native?unkonwn?) ****")}
  if (verbose & excludeUnknownRanges) {print('INFO: parameters set so records in unknown ranges are filtered here. Make sure this is what you want')}
  if (verbose & excludeNotmatchCountry) {print('INFO: parameters set so records that do not match recorded country vs. coordinate countries are filtered here Make sure this is what you want')}
  Analysis.F <-occProfileR::countryStatusRangeAnalysis(df=dat,
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
                                                       verbose = T)
  
  if  (nrow (Analysis.F$stay) != 0 ) {   dat.Q.F <- Analysis.F$stay }#; dat.Q.F$quality.grade <- 'F'
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
  tictoc:::toc()
  ############################################################################
  ### STEP 7:  Quality A-E Environmental and Geographical outliers  - analysis chunk
  ############################################################################
  #set timer
  tictoc:::tic('Total Geographic and Range Analysis:')
  print('Total Geographic and Range Analysis started....')
  
  if (verbose) {print ("**** RESOLIVNG  ****")}
  
  #ANALYSIS ELEMENTS
  #this is important for development, need to specify the number of ELEMENTS of analysis
  #to sumarize the results later on we will need that number
  N_analysis = 8
  
  ### ELEMENT 1: CENTROID ISSUE DETECTION
  tictoc:::tic('Centroid detection')
  print ('Centroid detection started ...')

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
                                       cfsf=countryfield.shapefile,
                                       method = methodCentroidDetection,
                                       do= doCentroidDetection)
  tictoc:::toc()
  
  ### ELEMENT 2: HYPER-HUMAN ENVIRONMENT
  tictoc:::tic('Land Use Land Cover analysis')
  print ('Land Use Land Cover analysis started ...')
  Analysis.2 <- humanDetection     (df = dat,
                                    xf = x.field,
                                    yf = y.field,
                                    .points.proj4string =points.proj4string,
                                    ras.hii = ras.hii,
                                    .th.human.influence =th.human.influence,
                                    do = doHumanDetection,output.dir=output.dir)
  tictoc:::toc()
  
  ### ELEMENT 3: BOTANICAL GARDEN PLACEMENT -- FROM LOCALITY NAME
  tictoc:::tic('Institution locality')
  print  ('Institution locality started ...')
  Analysis.3 <- institutionLocality (df=dat,xf = x.field,yf=y.field,
                                     lf=l.field,
                                     do = doInstitutionLocality,
                                     method = methodInstitutionLocality)
  tictoc:::toc()
  
  ### ELEMENT 4: GEOGRAPHICAL OUTLIER
  tictoc:::tic('geographic outliers detection')
  print ('geographic outliers detection started')
  Analysis.4 <- geoOutliers(df=dat,
                            xf=x.field,
                            yf=y.field,
                            .alpha.parameter = alpha.parameter,
                            #do=doGeoOutliers,
                            method = methodGeoOutliers,
                            .projString = points.proj4string)
  tictoc:::toc()
  
  ### ELEMENT 5: ADD the geoIndicator of wrong country if that has not been determined to be an exclusive dimension
  # tictoc:::tic('Country geoindication')
  # if (!excludeUnknownRanges & doRangeAnalysis){
  #   
  #   a = dat$countryStatusRange_wrongNTV_test ; a[is.na(a)] <- 0
  #   b = dat$countryStatusRange_wrongINV_test ; b[is.na(b)] <- 0
  #   
  #   Analysis.5 <-  a+b
  #   Analysis.5 <-  ifelse(test = Analysis.5>1,yes = 1,no = 0)
  #   Analysis.5 <- data.frame (countryStatusRange_notinKnownRange_test=Analysis.5, countryStatusRange_notinKnownRange_comments='')
  #   Analysis.5$unknownRange_score <- occProfileR:::.gimme.score(Analysis.5)
  #   row.names(Analysis.5) <- NULL
  #   rm(a);rm(b)
  # }
  # if (excludeUnknownRanges | !doRangeAnalysis){
  #   Analysis.5 <- data.frame (countryStatusRange_notinKnownRange=NA,
  #                             countryStatusRange_notinKnownRange_comments='unkonown ranges are classifed as F directly',
  #                             unknownRange_score=NA)[1:nrow(dat),]
  #   row.names(Analysis.5) <- NULL
  # }
  # 
  # ### ELEMENT 6: ADD the geoIndicator of country recorded different than country coordinates
  # if (!excludeNotmatchCountry & doCountryRecordAnalysis){
  #   Analysis.6 = dat$countryStatusRange_wrongCTRY_test
  #   Analysis.6 <- data.frame (countryStatusRange_wrongReportCtry_test=Analysis.6,
  #                             countryStatusRange_wrongReportCtry_comments='')
  #   Analysis.6$wrongReportCtry_score <- occProfileR:::.gimme.score(Analysis.6)
  # }
  # if (excludeNotmatchCountry | !doCountryRecordAnalysis){
  #   Analysis.6 <- data.frame (countryStatusRange_wrongReportCtry=NA,
  #                             countryStatusRange_wrongReportCtry_comments='Analysis not performed',
  #                             wrongReportCtry_score= NA)[1:nrow(dat),]
  #   row.names(Analysis.6) <- NULL
  # }
  # tictoc:::toc()
  
  ### ELEMENT 7: ENVIRONMENTAL OUTLIER
  tictoc:::tic('Environmental outliers')
  print ('Environmental outliers analysis started...')
  Analysis.7 <- envOutliers (.r.env=r.env,
                             df= dat, xf=x.field,
                             yf =y.field,
                             .th.perc.outenv = th.perc.outenv,
                             .sp.name = sp.name,
                             .projString = points.proj4string,
                             method = methodEnvOutliers,
                             do = doEnvOutliers)
  
  tictoc:::toc()
  
  ### ELEMENT 8: Coordinate accuracy
  tictoc:::tic('geoEnvironmental accuracy')
  print ('geoEnvironmental accuracy analysis started...')
  
  Analysis.8 <- geoEnvAccuracy(df=dat,
                               xf = x.field,
                               yf = y.field,
                               af = a.field,
                               dsf= ds.field,
                               r.env = r.env,
                               ef= e.field,
                               raster.elevation = r.dem,
                               do = do.geoEnvAccuracy,method = methodGeoEnvAccuracy,doParallel=doParallel,mc.cores=mc.cores)
  tictoc:::toc()
  
  ### SUMMARY ANALYSIS RESULTS
  list.analysis = list()
  for (i in 1:N_analysis){
    if (exists(paste0('Analysis.',i))){
      list.analysis[[i]] <- get (paste0('Analysis.',i))
    } else {list.analysis[[i]] =  NULL}
  }
  
  df.qualityAssessment <- dplyr::bind_cols(list.analysis)
  row.names(df.qualityAssessment) <- NULL
  
  
  #timer for the analytic processes
  tictoc:::toc()
  
  ############################################################################
  #### STEP 8:  Quality A-E Environmental and Geographical outliers
  # assignation of grade to a record
  ###########################################################################
  # tictoc:::tic('Labeling occurrences')
  # 
  # if (gradingSettings$grading.test.type == 'strict')   {test.strictness.value = 0 }
  # if (gradingSettings$grading.test.type == 'majority') {test.strictness.value = 0.5 }
  # if (gradingSettings$grading.test.type == 'relaxed')  {test.strictness.value = 0.6 }
  # #start assigning values depending on conditions
  # suppressMessages (attach(df.qualityAssessment))
  # df.qualityAssessment$quality.grade <- NA
  # 
  # #load lazylogic function
  # .lazylogic <- function (e,...){
  #   
  #   o <- eval(parse(text=e))
  #   ifelse(is.na(o),F,o)
  #   
  # }
  # 
  # 
  # #E is set so that it could be optional
  # if ('E' %in% gradingSettings$qualifier.label.scoping){
  #   grade.E <- (.lazylogic (e = 'HumanDetection_score >= test.strictness.value')) | (.lazylogic (e = 'institutionLocality_score >= test.strictness.value')) | (.lazylogic (e = 'centroidDetection_score >= test.strictness.value')) | (.lazylogic (e = 'unknownRange_score >= test.strictness.value')) | (.lazylogic (e = 'wrongReportCtry_score >= test.strictness.value') | (.lazylogic (e = 'envOutliers_missingEnv_score >= test.strictness.value')))
  #   df.qualityAssessment$quality.grade[grade.E] <- 'E'
  # }
  # to.continue <- is.na(df.qualityAssessment$quality.grade)
  # 
  # grade.D <- (.lazylogic (e = 'geoOutliers_score >= test.strictness.value') ) &  (.lazylogic ('envOutliers_score  >= test.strictness.value') )
  # df.qualityAssessment$quality.grade[to.continue & grade.D] <- 'D'
  # to.continue <- is.na(df.qualityAssessment$quality.grade)
  # 
  # grade.C <- .lazylogic ('envOutliers_score  >= test.strictness.value')
  # df.qualityAssessment$quality.grade[to.continue & grade.C] <- 'C'
  # to.continue <- is.na(df.qualityAssessment$quality.grade)
  # 
  # grade.B <- .lazylogic (e = 'geoOutliers_score >= test.strictness.value') | (.lazylogic (e = 'geoenvLowAccuracy_score >= test.strictness.value') ) 
  # df.qualityAssessment$quality.grade[to.continue & grade.B] <- 'B'
  # to.continue <- is.na(df.qualityAssessment$quality.grade)
  # 
  # df.qualityAssessment$quality.grade[to.continue] <- 'A'
  # 
  # detach(df.qualityAssessment)
  # 
  # #add labels to the data
  # dat <- cbind(dat,df.qualityAssessment)
  # tictoc:::toc()
  # 
  # ############################################################################
  ### STEP 9: BUILD FULL dataframe
  ############################################################################
  previous.Q.objects <- grep(pattern = 'dat.Q.',ls(),value = T)
  full.qaqc <- cbind(dat, df.qualityAssessment)
  
  for (o in previous.Q.objects){
    rowsToAdd = get(o)
    if (nrow (rowsToAdd)> 0) {full.qaqc <- dplyr::bind_rows(full.qaqc, get(o))}
    
  }
  
  ###########################################################################
  ### STEP 10: GET QUALIFIERS (attributes you want to add to your label that help decide its use)
  ###########################################################################
  # tictoc:::tic('Adding qualifiers name tags to  occurrences')
  # 
  # #check if user does not want the qualifiers, and return output
  # if (gradingSettings$qualifiers == F) {
  #   
  #   full.qaqc$qualifiers <- NA
  #   full.qaqc$quality.label <- full.qaqc$quality.grade
  #   
  #   #WRITE OUTPUTS AND END THE PROCESS
  #   if (write.full.output==T) {
  #     write.csv (full.qaqc,   paste0(output.dir,'/',sp,'_long_' ,
  #                                    output.base.filename,'.csv'),row.names = F)
  #   }
  #   
  #   short.qaqc<-full.qaqc[,c(x.field,y.field,'quality.grade',
  #                            'qualifiers','quality.label')]
  #   if (write.simple.output==T) {
  #     write.csv (short.qaqc,
  #                paste0(output.dir,'/',sp,'_short_' ,output.base.filename,
  #                       '.csv'),row.names = F)
  #   }
  #   
  #   output.function <- list (occ_full_profile=full.qaqc,
  #                            occ_short_profile=short.qaqc)
  #   return (output.function)
  # }
  # 
  # #get the qualifiers
  # tag1 <- occProfileR:::.qualifier.invasive.range(df=full.qaqc,
  #                                                 tag = 'i',
  #                                                 qualifier.lab.scp = qualifier.label.scoping)
  # tag2 <- occProfileR:::.qualifier.native.range(df=full.qaqc,
  #                                               tag = 'n',
  #                                               qualifier.lab.scp = qualifier.label.scoping)
  # 
  # #we will also add if we have a time stamp
  # tag3 <-  occProfileR:::.qualifier.timestamp(df=full.qaqc,
  #                                             tf = t.field,
  #                                             tag = 'T',
  #                                             qualifier.lab.scp = qualifier.label.scoping)
  # 
  # #we will also add if we have the right precision
  # tag4 <- occProfileR:::.qualifier.precision(df=full.qaqc,tag='p',
  #                                            .s = res.in.minutes,
  #                                            .e = res.in.minutes,
  #                                            xf = x.field,
  #                                            yf = y.field,
  #                                            qualifier.lab.scp = qualifier.label.scoping)
  # 
  # #we will also add an elevation threshold level
  # tag5 <- occProfileR:::.qualifier.elevation (df=full.qaqc,tag='e',
  #                                             qualifier.lab.scp = qualifier.label.scoping 
  # )
  # 
  # #get all tags together
  # tags.applied <- grep('tag',ls(),value=T)
  # tag.output <- lapply (1:length(tag1), function (i){
  #   
  #   row.tag <- unlist(lapply (tags.applied, function (x){a <- get(x)[i] }))
  #   row.tag <- paste(row.tag,collapse = '-')
  #   
  #   return(row.tag)
  #   
  # })
  # tag.output <- unlist (tag.output)
  # tag.output <- lapply (tag.output, function (x){
  #   
  #   row.tag <- strsplit (x,split='-')[[1]]
  #   row.tag <- row.tag [row.tag!="NA"]
  #   if (length(row.tag)==0){return(NA)}
  #   if (length(row.tag)==1){return(row.tag)}
  #   if (length(row.tag)>1) {row.tag <- paste(row.tag,collapse = '-'); return (row.tag)}
  #   
  # })
  # tag.output <- unlist (tag.output)
  # full.qaqc$qualifiers <- tag.output
  # 
  # #collapse to a unique label
  # full.qaqc$quality.label <- occProfileR:::.paste3 (as.character(full.qaqc$quality.grade),full.qaqc$qualifiers,sep = '/')
  # 
  # #reorder as the same as inputs
  # full.qaqc <- full.qaqc [sort (full.qaqc$roworder,decreasing = F, na.last=T),]
  # 
  # tictoc:::toc()
  # 
  ########################################################################
  ### STEP 11: WRITE THE OUTPUTS
  ########################################################################
  tictoc:::tic('Preparing and Writing outputs')
  print('Preparing and Writing outputs started ...')
  #reorder data as original
  full.qaqc = full.qaqc[order (full.qaqc$roworder),]
  full.qaqc= full.qaqc[,! names (full.qaqc)== 'roworder']
  
  #write outputs
  if (write.full.output==T) {
    sp2 = occProfileR:::.join.spname (sp)
    newdir = paste0(output.dir,'/',sp2)
    dir.create (newdir,recursive = T,showWarnings = F)
    written = try (write.csv (full.qaqc,  paste0(newdir,'/',
                                                 output.base.filename,'_',sp,'_long.csv'),
                              row.names = F),silent = T)
    if (class(written)=='try-error') save (list = 'full.qaqc',file = paste0(newdir,'/',output.base.filename,'_',sp,'_long.RData'))
    if (class(written)=='try-error') try (file.remove(paste0(newdir,'/',output.base.filename,'_',sp,'_long.csv')), silent=T )
  }
  
  #short qaqc
  idCols1 = which (names (full.qaqc)== x.field )
  idCols2 = which (names (full.qaqc)== y.field )
  idCols3 = grep (pattern = '_test',names(full.qaqc))
  idCols4 = grep (pattern = '_score',names(full.qaqc))
  short.qaqc = full.qaqc[,c(idCols1,idCols2,idCols3,idCols4)]
  
  if (write.simple.output==T) {
    sp2 = occProfileR:::.join.spname (sp)
    newdir = paste0(output.dir,'/',sp2)
    dir.create (newdir,recursive = T,showWarnings = F)
    write.csv (short.qaqc,
               paste0(newdir,'/',output.base.filename,'_',sp,
                      '_short.csv'),row.names = F)
  }
  tictoc:::toc()
  
  output.function <- list (occTest_full=full.qaqc, occTest_short=short.qaqc)
  tictoc:::toc()
  
  return (output.function)
  
  
}

