# occTest ========================================

#'Occurrence tests 
#'
#'Perform tests for data quality in species occurrence  using several methods
#' @param sp.name character. Name of the species.
#' @param habitat NULL
#' @param sp.table data.frame. Object with the coordinate data.
#' @param r.env raster or rasterStack. Environmental data(e.g. typically climatic.
#' @param tableSettings list. Elements corresponding to different settings of the input occurrence table. 
#' @param analysisSettings list. Elements corresponding to different settings of the analysis functions . 
#' @param gradingSettings list. Not implemented yet. Defaults to NULL.Elements corresponding to different settings of the analysis functions . 
#' @param writeoutSettings list. Elements corresponding to different settings of the analysis functions . 
#' @param return.spatial.data logical. Should the spatial dataset of \code{analysisSettings} attached to the metadata?, default is FALSE to save memory
#' @param r.dem raster. Elevation data (in meters).
#' @param ntv.ctry character. vector with ISO3 code of the countries where species is considered native
#' @param inv.ctry character. vector with ISO3 code of the countries where species is considered invasive 
#' @param resolveAlienCtry logical. To automatically try to detect  countries for which species is considered native
#' @param resolveNativeCtry logical. To automatically try to detect  countries for which species is considered alien
#' @param interactiveMode logical. Should prompts be ouput for some decisions taken by the workflow? Deafult FALSE,
#' @param verbose logical. Print workflow information? Default to FALSE
#' @param doParallel logical. Should some operations be run in parellel when possible? Default to FALSE
#' @param mc.cores numeric. If doParallel is TRUE, then how many cores to use? Default to 2
#' @return data frame with the tests performed (field $_test), specific comment for the tests ($_comments), the exact value of the test ($_value), and scores summarizing results across tests with an objective ($_score) 
#' @note 
#' The output dataframe allows users to classify or scrub the occurrences the way they want to. \cr
#' The names of the columns in the output object have the following naming convention: \cr
#' $AnalysisType$_$SpecificTest$_value: numeric or logical. Shows the quantitative result of the test (sometimes the same as in the result of the _test) \cr
#' $AnalysisType$_$SpecificTest$_test: logical Shows whether the occurrence passes or not the test, being TRUE a flag for a wrong record and NA indicating that the test was not performed on that record. \cr
#' $AnalysisType$_$SpecificTest$_comment: character. Shows some comments related to the specific test.\cr
#' Examples: HumanDetection_HumanInfluence_value gives you the score of current human influence in the record HumanDetection_HumanInfluence_test gives you whether we consider the former value an error/bias (TRUE) or not (FALSE) HumanDetection_HumanInfluence_comment gives you a commen that give further detail on the analysis. In this case that the threshold of 45 was used for the test. HumanDetection_score summarizes all the other HumanDetection tests and outputs a value from 0 to 1. A value of 0.5 would indicate that half of the tests used indicate that is an a Human signal in the record.
#' @examples \donttest{
#' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#' #load environmental raster
#' library (raster)
#' library (sf)
#' library (occTest)
#' #load occurrence data
#' occData <- read.csv (system.file('ext/exampleOccData.csv',package = 'occTest'))
#' #load environmental raster
#' renv <- raster (system.file('ext/AllEnv.tif',package = 'occTest'))
#' #load elevation raster
#' dem <- raster (system.file('ext/DEM.tif',package = 'occTest'))
#' #load settings
#' settings <- readRDS (system.file('ext/exSettings.rds',package = 'occTest'))
#' #run occTest
#' out = occTest(sp.name='MyFake species',
#'              sp.table = occData,ntv.ctry = 'ESP',inv.ctry = 'FRA',
#'              tableSettings = settings$tableSettings,
#'              writeoutSettings = settings$writeoutSettings,
#'              analysisSettings = settings$analysisSettings,
#'              r.env = renv,r.dem=dem)
#' }
#' @export
occTest = function(
  sp.name,
  habitat=NULL,
  sp.table,
  r.env,
  
  tableSettings=NULL,
  analysisSettings=NULL,
  writeoutSettings=NULL,
  gradingSettings=NULL,
  return.spatial.data=FALSE,
  
  r.dem=NULL,
  ntv.ctry=NULL,
  inv.ctry=NULL,
  resolveAlienCtry=FALSE,
  resolveNativeCtry=FALSE,
  
  interactiveMode=FALSE,
  verbose = FALSE,
  doParallel=FALSE,
  mc.cores=2){

  tictoc::tic()
  
  ### STEP 00: Initial checks ====
  #set timer
  tictoc::tic('Initial checks and formatting')
  message('Initial checks and formatting started...')
  
  #identify starting issues and convert to the right type of object
  if(missing(sp.table)) stop('missing sp.table')
  if(missing(r.env)) stop('missing r.env')
  if(! pingr::is_online()) stop('You seem not to have Internet connection. This package requires internet connection for several tests. Please go online')
  if(inherits(r.env,'SpatRaster')) stop ('Sorry, occTest is not ready for terra pkg yet. Transform your environmental raster to a rasterLayer')
  
  sp.table <- as.data.frame (sp.table)
    
  
  
  ### STEP 0: Load settings and study native and invasive countries ====
  if (is.null(tableSettings) | is.null(analysisSettings) | is.null(writeoutSettings)) defaultSettings = occTest::minimalSettings()
  
  #load table settings(old stuff, we  could attach the different labels)
  if(is.null(tableSettings)){ tableSettings = defaultSettings$tableSettings}
  
  x.field             = tableSettings$x.field
  y.field             = tableSettings$y.field
  if (! all(c(x.field,y.field)%in%  names(sp.table))){stop('No coordinate fields specified')}
  
  #set table field names
  ## taxonobservation.id
  taxonobservation.id = tableSettings$taxonobservation.id
  if (!is.null(taxonobservation.id)){
    if(! taxonobservation.id %in% names(sp.table))   taxonobservation.id = NULL
  }
  ## t.field
  t.field             = tableSettings$t.field
  if (!is.null(t.field)){
    if(! t.field %in% names(sp.table))   t.field = NULL
  }
  ## l.field
  l.field             = tableSettings$l.field  
  if (!is.null(l.field)){
    if(! l.field %in% names(sp.table))   l.field = NULL
  }
  ## c.field
  c.field             = tableSettings$c.field
  if (!is.null(c.field)){
    if(! c.field %in% names(sp.table))   c.field = NULL
  }
  ## e.field
  e.field             = tableSettings$e.field
  if (!is.null(e.field)){
    if(! e.field %in% names(sp.table))   e.field = NULL
  }
  ## a.field
  a.field             = tableSettings$a.field
  if (!is.null(a.field)){
    if(any(! a.field %in% names(sp.table)))   a.field = NULL
  }
  ## ds.field
  ds.field            =  tableSettings$ds.field
  if (!is.null(ds.field)){
    if(! ds.field %in% names(sp.table))   ds.field = NULL
  }

  
  #load analysisSettings
  if(is.null(analysisSettings)) analysisSettings = defaultSettings$analysisSettings
  coordinate.decimal.precision = analysisSettings$geoSettings$coordinate.decimal.precision
  points.proj4string = analysisSettings$geoSettings$points.proj4string
  doCoastalReassignment = analysisSettings$doCoastalReassignment
  countries.shapefile = analysisSettings$countryStatusRange$countries.shapefile
  landSurfacePol = analysisSettings$landSurfacePol
  if (any (class(countries.shapefile) %in% c('sf'))) {stop('sf not implemented, tramsform ctry polygons to Spatial object and rerun')}
  countryfield.shapefile = analysisSettings$countryStatusRange$countryfield.shapefile
  doRangeAnalysis = analysisSettings$countryStatusRange$doRangeAnalysis
  excludeUnknownRanges = analysisSettings$countryStatusRange$excludeUnknownRanges
  excludeNotmatchCountry = analysisSettings$countryStatusRange$excludeNotmatchCountry
  doCountryRecordAnalysis = analysisSettings$countryStatusRange$doCountryRecordAnalysis
  
  doCentroidDetection = analysisSettings$centroidDetection$doCentroidDetection
  methodCentroidDetection = analysisSettings$centroidDetection$methodCentroidDetection
  
  doHumanDetection = analysisSettings$humanDetection$doHumanDetection
  methodHumanDetection = analysisSettings$humanDetection$methodHumanDetection
  th.human.influence = analysisSettings$humanDetection$th.human.influence
  ras.hii = analysisSettings$humanDetection$ras.hii
  
  doLandUseSelect = analysisSettings$landUseType$doLandUse
  methodLandUseSelect = analysisSettings$landUseType$methodLandUse
  landUseCodes = analysisSettings$landUseType$landUseCodes
  ras.landUse = analysisSettings$landUseType$ras.landUse

  doInstitutionLocality = analysisSettings$institutionLocality$doInstitutionLocality
  methodInstitutionLocality = analysisSettings$institutionLocality$methodInstitutionLocality
  
  doGeoOutliers = analysisSettings$geoOutliers$doGeoOutliers
  methodGeoOutliers = analysisSettings$geoOutliers$methodGeoOutliers
  alpha.parameter = analysisSettings$geoOutliers$alpha.parameter
  
  doEnvOutliers = analysisSettings$envOutliers$doEnvOutliers
  methodEnvOutliers = analysisSettings$envOutliers$methodEnvOutliers
  th.perc.outenv= analysisSettings$envOutliers$th.perc.outenv
  
  methodGeoEnvAccuracy=analysisSettings$geoenvLowAccuracy$methodGeoEnvAccuracy
  doGeoEnvAccuracy=analysisSettings$geoenvLowAccuracy$doGeoEnvAccuracy
  elev.quality.threshold=analysisSettings$geoenvLowAccuracy$elev.quality.threshold
  
  #load gradingSettings
  # if(is.null(gradingSettings)){ gradingSettings = defaultSettings$gradingSettings}
  # 
  # grading.test.type = gradingSettings$grading.test.type
  # qualifiers = gradingSettings$qualifiers
  # qualifier.label.scoping = gradingSettings$qualifier.label.scoping
  
  #load writeOutSettings
  if(is.null(writeoutSettings)) writeoutSettings = defaultSettings$writeoutSettings
  output.dir = writeoutSettings$output.dir
  if(is.null(output.dir)) output.dir = tempdir(); message('Output directory not specified. Set to a temporary directory')
  writeAllOutput = writeoutSettings$writeAllOutput
  write.simple.output =  writeoutSettings$write.simple.output
  write.full.output = writeoutSettings$write.full.output
  output.base.filename = writeoutSettings$output.base.filename
  
  ### STEP 1a: Data formatting and compatibility for biogeo and initial checks =====

  #add fields necesary for initial table
  sp =  .join.spname(sp.name)
  sp.table$Species = sp
  
  sp.table2 =  .checkfields(dat=sp.table,xf = x.field, yf=y.field,
                                     ef = e.field,tf = t.field,lf = l.field,
                                     cf = c.field,          
                                     idf = taxonobservation.id)
  
  dat =  .addmainfields2(sp.table2,species = 'Species')
  dat$comments = rep('', nrow(dat))
  
  #check data structure
  ck  =  .checkdatastr2(dat,xf = x.field,yf=y.field)
  if(sum(ck$Present)!=10){stop("Error: required table fields could not be created")}
  
  #For development:
  #placeholder to convert other kind of coordinate data to numerical latlong
  #or whatever decided projection we want to work on
  
  #CHECK 4
  #ensure that all the geospatial data inputs are is in the same coordinate system
  potential.geosp.objects = list(countries.shapefile,
                                  r.env,
                                  r.dem,
                                  ras.hii,
                                  points.proj4string)
  x = sapply(potential.geosp.objects, is.null)
  actual.input.geosp.objects = potential.geosp.objects[-x]
  actual.input.geosp.objects =  .subsetlist.nonNULL(actual.input.geosp.objects)
  
  #the development should be in the direction of automatically check and transform
  # .check.geospatial.data(list.geospatial.objects =actual.input.geosp.objects)
  #lapply(actual.input.geosp.objects, function(x)raster::projection(x))
  
  
  #CHECK 5
  #get into precision of decimals
  for(i in c(x.field, y.field)){
    dat[,i] = round(dat[,i], digits=coordinate.decimal.precision)
  }
  
  #CHECK 6
  #need tweaking: this will only happen if the data is in geographic coordinates, because the unites will be in minutes
  #res is is in minutes...therefore(30 arcsec = 0.5 min)
  #if this is a projected data, we need to redo stuff and also tweak funcitons in biogeo
  #maybe add a warning or something in check geospatial data? idk
  res.in.minutes = raster::res(r.env[[1]])[1] * 60
  
  tictoc::toc()
  
  ### STEP 1b [OPTIONAL]: Automatically solve native or invasive range  =====
  #### NOT IMPLEMENTED YET !!!!!!!!! TO DEVELOP WITH BIEN GNRS/NRS ,...
  #whatever the want we are going to set it to False from now
  if (any (c(interactiveMode,resolveNativeCtry,resolveAlienCtry))) message('Automatic resolving of native and non-native ranges not implemented yet')
  interactiveMode= FALSE
  resolveNativeCtry = FALSE
  resolveAlienCtry = FALSE
  #set timer
  #tictoc::tic('Resolve native and invasive countries')
  
  #automatically resolve invasive and native countries for target species(not implemented yet)
  
  if(interactiveMode & is.null(ntv.ctry)){
    resolveNativeCtry = if(interactive()) utils::askYesNo(default = FALSE,msg = "You have not provided countries for the species native range. Do you want to infer from global databases?")
  }
  if(interactiveMode & is.null(inv.ctry)){
    resolveAlienCtry = if(interactive()) utils::askYesNo(default = FALSE,msg = "You have not provided countries for the species alien range. Do you want to infer from global databases?")
  }
  if(any(resolveNativeCtry,resolveAlienCtry)){
    if(verbose){message("**** RESOLIVNG NATIVE AND INVASIVE RANGES ****")}
    
    xydatTemporary = dat[,c(x.field,y.field)]
    xydatTemporary = xydatTemporary[stats::complete.cases(xydatTemporary),]
    
    
    checkCountries = occTest::nativeStatusCtry(spName = sp.name, xydat=xydatTemporary, resolveNative = resolveNativeCtry, resolveAlien = resolveAlienCtry ,verbose = verbose)
    if(resolveNativeCtry){
      ntv.ctry = c(ntv.ctry,checkCountries$ntvCtry)
      ntv.ctry = unique(ntv.ctry)
    }
    if(resolveAlienCtry){
      inv.ctry = c(inv.ctry,checkCountries$invCtry)
      inv.ctry = unique(inv.ctry)
    }
    
  }
  #tictoc::toc()
  
  if(verbose){message("**** FILTERING PHASE STARTED  ****")}
  
  ### STEP 2: Quality Filter : Identify records wo spatial info =====
  #set timer
  tictoc::tic('Filter major coordinate Issues')
  message('Filter major coordinate Issues started...')

  Analysis.H = filterMissing(df = dat,xf = x.field,yf = y.field)
  dat.Q.H = Analysis.H$stay
  dat = Analysis.H$continue
  
  #valid coordinates in geographic projections and zero zero issues anc decimal conversion issues
  if(as.character(points.proj4string)%in% 
      c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',"+proj=longlat +datum=WGS84 +no_defs")){
    
  
    coordIssues_invalidCoord_value =! CoordinateCleaner::cc_val(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=FALSE)
    dat$coordIssues_invalidCoord_value = coordIssues_invalidCoord_value
    coordIssues_invalidCoord_test = as.logical(coordIssues_invalidCoord_value)
    dat$coordIssues_invalidCoord_test = coordIssues_invalidCoord_test
    
    dat.Q.H1 = dat[coordIssues_invalidCoord_test,]
    dat = dat[!coordIssues_invalidCoord_test,]
    
    
    status.out= .status.tracker.and.escaping(dataset.to.continue = dat,
                                                      wfo = write.full.output,
                                                      wso = write.simple.output,
                                                      xf = x.field,
                                                      yf=y.field,
                                                      od = output.dir,
                                                      rsd=return.spatial.data,
                                                      obf = output.base.filename,
                                                      sp=sp,
                                                      as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
    if(any(class(status.out)=='occTest')) {return(status.out)}
    
    #zero zero long lat
    coordIssues_zeroCoord_value = ! CoordinateCleaner::cc_zero(x = dat,lon = x.field,lat = y.field,value='flagged',verbose=FALSE)
    dat$coordIssues_zeroCoord_value = coordIssues_zeroCoord_value
    coordIssues_zeroCoord_test = as.logical(coordIssues_zeroCoord_value)
    dat$coordIssues_zeroCoord_test = coordIssues_zeroCoord_test
    dat.Q.H2 = dat[coordIssues_zeroCoord_test,]
    dat = dat[!coordIssues_zeroCoord_test,]
    
    status.out= .status.tracker.and.escaping(dataset.to.continue = dat,
                                                      wfo = write.full.output,
                                                      wso = write.simple.output,
                                                      xf = x.field,
                                                      yf=y.field,
                                                      od = output.dir,
                                                      rsd=return.spatial.data,
                                                      obf = output.base.filename,
                                                      sp=sp,
                                                      as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
    if(any(class(status.out)=='occTest')) {return(status.out)}
    
    #error conversion of decimals 
    dat.tmp = dat 

    if(is.null(ds.field)){ warning('No dataset field provided, considering everything as a unique dataset')
      dat.tmp$MyInventedCommonDataset = 'TemporaryDatasetName' 
      ds.field = 'MyInventedCommonDataset'}
    outDecimalTest = try({ .cd_ddmm_occTest(x = dat.tmp,lon = x.field,lat = y.field,ds=ds.field , value='flagged',verbose=FALSE,diff=1.5)},silent=TRUE)
    if (class(outDecimalTest) %in% c('error','try-error'))     outDecimalTest = rep (NA,length.out=nrow (dat.tmp))
    coordIssues_coordConv_value = ! outDecimalTest
    dat$coordIssues_coordConv_value = coordIssues_coordConv_value
    coordIssues_coordConv_test = as.logical(coordIssues_coordConv_value)
    dat$coordIssues_coordConv_test = coordIssues_coordConv_test
    
    if (all(is.na(coordIssues_coordConv_test)))      warning ('Coordinate Conversion test errors not perfomred due to internal error. Continuing without that test')
    
    if (!all(is.na(coordIssues_coordConv_test)))     dat.Q.H3 = dat[coordIssues_coordConv_test,]
    if (!all(is.na(coordIssues_coordConv_test)))     dat = dat[!coordIssues_coordConv_test,]

    if(ds.field == 'MyInventedCommonDataset')ds.field = NULL
    rm(dat.tmp)
    
    status.out= .status.tracker.and.escaping(dataset.to.continue = dat,
                                                      wfo = write.full.output,
                                                      wso = write.simple.output,
                                                      xf = x.field,
                                                      yf=y.field,
                                                      od = output.dir,
                                                      rsd=return.spatial.data,
                                                      obf = output.base.filename,
                                                      sp=sp,
                                                      as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
  }
  
  if(any(class(status.out)=='occTest')) {return(status.out)}
  
  #indicate issues of georeference and put them aside
  obj.issues = c('dat.Q.H','dat.Q.H1','dat.Q.H2','dat.Q.H3')
  obj.exist = sapply(obj.issues ,FUN = function(x)exists(x))
  obj.issues = obj.issues[obj.exist]
  if(any(obj.issues %in% c('dat.Q.H1','dat.Q.H2','dat.Q.H3'))){
    if(length(obj.issues)>0){
      dat.excl.H = lapply(obj.issues, function(x){
        
        if(nrow(get(x))> 0 ){
          a = get(x)
          #a$quality.grade = 'H'
          a$Exclude = 1 
          if(x == 'dat.Q.H')
            a$Reason = "No coordinates"
          if(x == 'dat.Q.H1')
            a$Reason = "No valid coords"
          if(x == 'dat.Q.H2')
            a$Reason = "LonLat 0 0"
          if(x == 'dat.Q.H3')
            a$Reason = "Coord decimal conversion error"
          return(a)
          
        } else { return(NULL)}
        
      })
      dat.excl.H = dplyr::bind_rows(dat.excl.H)
    }
    if(exists('dat.excl.H')){
      if(!is.null(get('dat.excl.H'))){
        if(nrow(dat.excl.H)>0){
          dat.Q.H = dat.excl.H
          rm(dat.excl.H)} else {rm(dat.excl.H)}
    }
    if(exists('dat.Q.H1')){rm(dat.Q.H1)} 
    if(exists('dat.Q.H2')){rm(dat.Q.H2)} 
    if(exists('dat.Q.H3')){rm(dat.Q.H3)} 
    
    }
    }
  
  #check outputs and escape ifneedbe //
  status.out= .status.tracker.and.escaping(dataset.to.continue = dat,
                                                    wfo = write.full.output,
                                                    wso = write.simple.output,
                                                    xf = x.field,
                                                    yf=y.field,
                                                    od = output.dir,
                                                    rsd=return.spatial.data,
                                                    obf = output.base.filename,
                                                    sp=sp,
                                                    as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
  
  if(any(class(status.out)=='occTest')) {return(status.out)}
  tictoc::toc()
  ### STEP 4: Filter Quality G : Identify duplicate records(in geographic space)to prevent pseudoreplicaton ===== 
  #set timer
  tictoc::tic('Filter duplicates')
  message('Filter duplicates started')
  if(verbose){message("**** RESOLVING: duplicates ****")}
  
  #indicate duplicates with exact coordinates
  Analysis.G = duplicatesexcludeAnalysis(df = dat,
                                         xf = x.field,
                                         yf = y.field,
                                         resolution.in.minutes=res.in.minutes,
                                         raster.grid = r.env[[1]])
  
  #here we merge together the two kind of duplicates, but maybe we would like to keep the relative duplicates for other purposes
  dat.Q.G = dplyr::bind_rows(Analysis.G$Dups.Grid, Analysis.G$Dups.Exact)
  if(nrow(dat.Q.G)== 0 ){rm(dat.Q.G)}
  #if (exists('dat.Q.G')){ dat.Q.G$quality.grade = 'G'}
  
  dat = Analysis.G$continue
  
  #check outputs and escape ifneedbe //
  status.out= .status.tracker.and.escaping(
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    rsd=return.spatial.data,
    obf = output.base.filename,
    sp=sp,
    as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
  
  
  if(any(class(status.out)=='occTest')) {return(status.out)}
  
  tictoc::toc()
  ### STEP 5: SEA/TERRESTRIAL POTENTIAL REASSIGNMENT AND RECHECK DUPLICATES  =====
  #-- and potential for new duplicates
  #set timer
  tictoc::tic('Resolving coastal reassignment')
  message('Resolving coastal reassignment started...')
  if(verbose){message("**** RESOLVING : sea/terrestrial reassignment ****")}
  #analysis of nearest cell next to the sea
  if (doCoastalReassignment) dat =  .nearestcell3(dat=dat,rst = r.env, xf=x.field, yf=y.field)
  
  #check results and recheck dups ifneedbe
  if(inherits(dat,'list')){
    dat = dat[[1]]
    moved.points = dat[['moved']]
    if(!is.null(output.dir)){
      sp.name2=  .join.spname(sp.name)
      odir = paste0(output.dir,'/',sp.name2)
      dir.create(odir,showWarnings = FALSE,recursive = TRUE)
      utils::write.csv( moved.points,
                 paste0(odir,'/',sp.name2,'_coastal_Reassignment.csv'))
    }
    
    
    
    ### RECHECK POTENTIAL DUPLICATES AGAIN AFTER REASSIGNATION
    Analysis.G.second.time=duplicatesexcludeAnalysis(df = dat,
                                                     xf = x.field,
                                                     yf = y.field,
                                                     resolution.in.minutes =
                                                       res.in.minutes)
    dat.Q.G.second.time = rbind(Analysis.G.second.time$Dups.Grid,
                                  Analysis.G.second.time$Dups.Exact)
    if (nrow(dat.Q.G.second.time)!=  0){
      #dat.Q.G.second.time$quality.grade = 'G'
      dat.Q.G = rbind(dat.Q.G,dat.Q.G.second.time)
      dat = Analysis.G.second.time$continue
    }
    
    rm(dat.Q.G.second.time)
    
  }
  ### filter those that are still in the sea/land (depending on habitatType)
  Analysis.LandSea = occTest::landSeaFilter(df = dat, 
                                            xf= x.field, y= y.field, 
                                            habType = habitat,
                                            habPol = landSurfacePol,
                                            verbose=verbose) 
  
  
  if (exists ('dat.Q.G') & nrow (Analysis.LandSea$stay) != 0)   dat.Q.G = dplyr::bind_rows (dat.Q.G, Analysis.LandSea$stay)
  if (!exists ('dat.Q.G') & nrow (Analysis.LandSea$stay) != 0)   dat.Q.G =  Analysis.LandSea$stay
  
  dat = Analysis.LandSea$continue
  
  #check outputs and escape if need be //
  status.out =  .status.tracker.and.escaping(
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    rsd=return.spatial.data,
    obf = output.base.filename,
    sp=sp,
    as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
  
  
  if(any(class(status.out)=='occTest')) {return(status.out)}
  tictoc::toc()
  
  ### STEP 6: Filter Quality : Country selection ======
  #set timer
  tictoc::tic('Resolving countryStatusRange Analysis')
  
  if(verbose) message('Resolving countryStatusRange Analysis started...')
  if(verbose & excludeUnknownRanges) message('INFO: parameters set so records in unknown ranges are filtered here. Make sure this is what you want')
  if(verbose & excludeNotmatchCountry) message('INFO: parameters set so records that do not match recorded country vs. coordinate countries are filtered here Make sure this is what you want')
  Analysis.F=countryStatusRangeAnalysis(df=dat,
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
                                        excludeNotmatchCountry =
                                          excludeNotmatchCountry,
                                        doRangeAnalysis = doRangeAnalysis,
                                        verbose = TRUE)
  
  if (nrow(Analysis.F$stay)!= 0 ) dat.Q.FALSE = Analysis.F$stay 
  dat = Analysis.F$continue
  
  #check outputs and escape ifneedbe //
  status.out= .status.tracker.and.escaping(
    dataset.to.continue = dat,
    wfo = write.full.output,
    wso = write.simple.output,
    xf = x.field,
    yf=y.field,
    od = output.dir,
    rsd=return.spatial.data,
    obf = output.base.filename,
    sp=sp,
    as=analysisSettings,ws = writeoutSettings,ts =tableSettings)
  
  
  if(any(class(status.out)=='occTest')) {return(status.out)}
  tictoc::toc()
  ### STEP 7: Environmental and Geographical outliers  - analysis chunk =====
  #set timer
  tictoc::tic('Analysis phase:')
  if(verbose){message("**** ANALYSIS PHASE STARTED  ****")}
  #ANALYSIS ELEMENTS
  ### ELEMENT : CENTROID ISSUE DETECTION
  tictoc::tic('Centroid detection')
  if(verbose) message('Centroid detection started ...')

  Analysis.1 = centroidDetection(.r.env = r.env,
                                 df = dat,
                                 xf = x.field,
                                 yf = y.field,
                                 cf = c.field,
                                 #to be changed according to definition of taxonobseration id
                                 #idf = taxonobservation.id,
                                 idf = NULL,
                                 .ntv.ctry = ntv.ctry,
                                 .inv.ctry = inv.ctry,
                                 .points.proj4string =
                                   points.proj4string,
                                 .countries.shapefile = countries.shapefile,
                                 cfsf=countryfield.shapefile,
                                 method = methodCentroidDetection,
                                 do= doCentroidDetection)
  tictoc::toc()
  
  ### ELEMENT : HYPER-HUMAN ENVIRONMENT
  tictoc::tic('Land Use Land Cover analysis')
  if(verbose) message('Land Use Land Cover analysis started ...')
  Analysis.2 = humanDetection (df = dat,
                               xf = x.field,
                               yf = y.field,
                               .points.proj4string =points.proj4string,
                               ras.hii = ras.hii,
                               .th.human.influence =th.human.influence,
                               do = doHumanDetection,output.dir=output.dir)
  tictoc::toc()
  
  ### ELEMENT : BOTANICAL GARDEN PLACEMENT -- FROM LOCALITY NAME
  tictoc::tic('Institution locality')
  if(verbose) message ('Institution locality started ...')
  Analysis.3 = institutionLocality(df=dat,xf = x.field,yf=y.field,
                                   lf=l.field,
                                   do = doInstitutionLocality,
                                   method = methodInstitutionLocality)
  
  tictoc::toc()
  
  ### ELEMENT : land use
  tictoc::tic('Records in land use')
  if(verbose) message ('Records in land use started...')
  Analysis.4 = landUseSelect(df=dat,xf = x.field,yf=y.field,
                             .points.proj4string =points.proj4string,
                             .landUseCodes = landUseCodes,
                             ras.landUse = ras.landUse,
                             method = methodLandUseSelect,
                             do = doLandUseSelect)
  
  tictoc::toc()
  
  ### ELEMENT : GEOGRAPHICAL OUTLIER
  tictoc::tic('geographic outliers detection')
  if(verbose) message('geographic outliers detection started')
  Analysis.5 = geoOutliers(df=dat,
                           xf=x.field,
                           yf=y.field,
                           .alpha.parameter = alpha.parameter,
                           do=doGeoOutliers,
                           method = methodGeoOutliers,
                           .projString = points.proj4string)
  tictoc::toc()
 
  ### ELEMENT 7: ENVIRONMENTAL OUTLIER
  tictoc::tic('Environmental outliers')
  if(verbose) message('Environmental outliers analysis started...')

  Analysis.6 = envOutliers(.r.env=r.env,
                             df= dat, xf=x.field,
                             yf =y.field,
                             .th.perc.outenv = th.perc.outenv,
                             .sp.name = sp.name,
                             .projString = points.proj4string,
                             method = methodEnvOutliers,
                             do = doEnvOutliers)
  tictoc::toc()
  
  ### ELEMENT 8: Coordinate accuracy
  tictoc::tic('geoEnvironmental accuracy')
  if(verbose) message('geoEnvironmental accuracy analysis started...')
  Analysis.7 = geoEnvAccuracy(df=dat,
                              xf = x.field,
                              yf = y.field,
                              af = a.field,
                              dsf= ds.field,
                              tf= t.field,
                              r.env = r.env,
                              ef= e.field,
                              raster.elevation = r.dem,
                              do = doGeoEnvAccuracy,
                              method = methodGeoEnvAccuracy,
                              doParallel=doParallel,
                              mc.cores=mc.cores)
  tictoc::toc()
  
  ### SUMMARY ANALYSIS RESULTS(NEED TO BE IMPROVED !)
  #this is important for development, need to specify the number of ELEMENTS of analysis
  #to sumarize the results later on we will need that number
  N_analysis = 8
  list.analysis = list()

    for(i in 1:N_analysis){
    if(exists(paste0('Analysis.',i))){
      list.analysis[[i]] = get(paste0('Analysis.',i))
    } else {list.analysis[[i]] =  NULL}
  }
  
  df.qualityAssessment = dplyr::bind_cols(list.analysis)
  row.names(df.qualityAssessment)= NULL
  
  #timer for the analytic processes
  tictoc::toc()

  ### STEP 9: BUILD FULL dataframe ====
  #load previous filtered objects
  previousFiltered = grep(pattern = 'dat.Q.',ls(),value = TRUE)
  full.qaqc = cbind(dat, df.qualityAssessment)
  for(o in previousFiltered){
    rowsToAdd = get(o)
    if(nrow(rowsToAdd)> 0){full.qaqc = dplyr::bind_rows(full.qaqc, get(o))}
    
  }
  
  ### STEP 10: WRITE THE OUTPUTS =====
  tictoc::tic('Preparing and Writing outputs')
  if(verbose) message('Preparing and Writing outputs started ...')
  #reorder data as original
  full.qaqc = full.qaqc[order(full.qaqc$roworder),]
  full.qaqc = full.qaqc[,! names(full.qaqc)== 'roworder']
  
  #write outputs
  if(write.full.output){
    sp2 =  .join.spname(sp)
    newdir = paste0(output.dir,'/',sp2)
    dir.create(newdir,recursive = TRUE,showWarnings = FALSE)
    written = try(utils::write.csv(full.qaqc,  
                            paste0(newdir,'/',output.base.filename,
                                   '_',sp,'_long.csv'),
                            row.names = FALSE),silent = TRUE)
    if(inherits(written,'try-error')) save(list = 'full.qaqc',file = paste0(newdir,'/',output.base.filename,'_',sp,'_long.RData'))
    if(inherits(written,'try-error')) try(file.remove(paste0(newdir,'/',output.base.filename,'_',sp,'_long.csv')), silent=TRUE )
  }
  tictoc::toc()
  
  #output.function = list(occTest_full=full.qaqc, occTest_short=short.qaqc)
  output.function = full.qaqc
  tictoc::toc()
  
  attr(output.function,"class")<-c("occTest",class(output.function))
 
  if(!return.spatial.data){
    analysisSettings$countryStatusRange$countries.shapefile<-NULL
    analysisSettings$humanDetection$ras.hii<-NULL
    analysisSettings$humanAnalysis$methodHyperHumanDetection<-NULL
    analysisSettings$rangeAnalysis$countries.shapefile<-NULL
  }
  
  
  attr(output.function,"Settings")<-list(tableSettings=tableSettings,analysisSettings=analysisSettings,writeoutSettings=writeoutSettings)

  
  return(output.function)
}

