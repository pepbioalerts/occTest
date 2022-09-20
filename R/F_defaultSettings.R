# defaultSettings ====
#' @title load default settings for occTest
#' @description Loads a list of lists with the different default parameters for analysis, outputs and grading needed in occTest
#' @details it can be use internally or it can be used by a user to subsequently modify parameters
#' @return list of lists with all different parameters to use in occProfile function
#' @keywords user
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
defaultSettings <- function (){
  
  #download info and files
  dest_url_hii = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/hii.rds'
  outFile_hii = paste0(tempdir(),'/hii.rds')
  if (!file.exists(outFile_hii)) utils::download.file(url=dest_url_hii,destfile = outFile_hii)
  
  
  
    defaultSettings = list (
    #grading Settings
    # gradingSettings = list (grading.test.type = 'majority', #other options are 'strict' 'relaxed'
    #                         qualifiers=T,
    #                         qualifier.label.scoping=c('A','B','C','D','E'))
    # 
    # ,
    #writing outputs settings
    writeoutSettings = list (#writing outputs
      output.dir=NULL,
      writeAllOutput=F, #overwrites write.simple.output, write.full.output
      write.simple.output=F,
      write.full.output=F,
      output.base.filename="occTest")
    ,
    #tableSettings
    tableSettings = list (taxonobservation.id = 'taxonobservationID',
                          x.field = 'decimalLongitude',
                          y.field = 'decimalLatitude',
                          t.field = 'eventDate', #time field (date)  
                          l.field = 'verbatimLocality', #locality field 
                          c.field = 'countryCode', #country field field (date)   
                          e.field = 'recordedElevation', #elevation recoreded  in meters
                          a.field = 'coordinateUncertaintyInMeters', #coordinate uncertainty in meters
                          ds.field = 'datasetName' #dataset field identifier
    )
    ,                     
    #analysis settings
    analysisSettings =list (
      geoSettings = list (
        coordinate.decimal.precision = 4,
        points.proj4string=sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      )
      ,
      countryStatusRange = list (
        countries.shapefile=rnaturalearth::ne_countries(scale = 50),
        countryfield.shapefile = 'iso_a3',
        doRangeAnalysis=T,
        excludeUnknownRanges= F,
        excludeNotmatchCountry= F,
        doCountryRecordAnalysis=T
      )
      ,
      centroidDetection = list (doCentroidDetection=T,
                               methodCentroidDetection='all'
      )
      ,
      humanDetection= list (doHumanDetection=T,
                           methodHumanDetection='all',
                           th.human.influence = 45,
                           ras.hii=readRDS(outFile_hii)

      )
      ,
      landUseType   = list (doLandUse=T,
                             methodLandUse='in',
                             landUseCodes = NULL,
                             ras.landUse=NULL #we need a default here to be downloaded
      )
      ,
      institutionLocality = list (doInstitutionLocality=T,
                                  methodInstitutionLocality='all'
      )
      ,
      
      geoOutliers = list (
        doGeoOutliers=T,
        methodGeoOutliers='all',
        alpha.parameter = 2
      )
      ,
      
      envOutliers = list (
        doEnvOutliers=T,
        methodEnvOutliers='all',
        th.perc.outenv =  0.2
      )
      ,
      geoenvLowAccuracy = list (methodGeoEnvAccuracy='all',
                               doGeoEnvAccuracy=T,
                               elev.quality.threshold = 100
      )
      
      
      
    )#end analysis settings
    
  )#end all default settings
  
  
  return (defaultSettings)
  
  
  
}

# showTableNames ====
#' @title Print naming conventions in occTest
#' @description prints a table with the the conventions used for column names
#' @details The function prints a guide to column naming conventions used by occTest in their default parameters. These defaults can be changed via setTableNames, but the user may also decide to format their input table according to these naming conventions. 
#' @return prints a dataframe
#' @keywords user
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
showTableNames <- function (){
  tabNames=readRDS(system.file('ext/tableColumns.rds',package='occTest'))
  DT::datatable(tabNames)
}

# setTableNames ====
#' @title set table names internally
#' @param x.field character. Name of the x coordinate field.
#' @param y.field character. Name of the y coordinate field.
#' @param t.field character. Name of the timestamp field.
#' @param l.field character. Name of the locality field.
#' @param c.field character. Name of the country code field.
#' @param e.field character. Name of the eleveation field.
#' @param a.field character. Name of the accuracy field.
#' @param ds.field character. Name of the dataset identifier field.
#' @param taxonobservation.id character. Name of the taxon observartion id field. 
#' @description
#' @details 
#' @return list
#' @keywords user
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
setTableNames <- function (x.field = NULL,
                           y.field = NULL,
                           t.field = NULL,
                           l.field = NULL,
                           c.field = NULL,
                           e.field = NULL,
                           a.field = NULL,
                           ds.field= NULL,
                           taxonobservation.id= NULL){
  targetList = defaultSettings()
  targetList = targetList$tableSettings
  ids = names (targetList)
  for (i in ids){
    if (!is.null(get(i))) targetList[[i]] <- get (i)
  }
  
  return (targetList)
}

# setTestTypes ====
#' @title Set the tests to run
#' @description function used to select which types of tests you want in occTest workflow
#' @details See occTest::showTests for further information on tests used in the packages
#' @param countryStatusRange logical. Should this test type be performed?
#' @param centroidDetection logical. Should this test type be performed?
#' @param humanDetection logical. Should this test type be performed?
#' @param landUseType logical. Should this test type be performed?
#' @param institutionLocality logical. Should this test type be performed?
#' @param geoOutliers logical. Should this test type be performed?
#' @param envOutliers logical. Should this test type be performed?
#' @param geoenvLowAccuracy logical. Should this test type be performed?
#' @return list with user modified settings
#' @keywords user
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
setTestTypes <- function (countryStatusRange = T,
                      centroidDetection = T,
                      humanDetection = T,
                      landUseType = T,
                      institutionLocality =T,
                      geoOutliers = T,
                      envOutliers = T,
                      geoenvLowAccuracy = T) {
  
  myTestDoParams= ls()
  targetList = defaultSettings()
  targetList = targetList$analysisSettings
  

  for (i in myTestDoParams){
    a = targetList[[i]]
    a[[grep(names(a),pattern = '^do')]] <-  get(i)
    targetList[[i]] = a
  }
  
  return (targetList)
}

# setTestBlocks ====
#' @title Set the tests to run
#' @description function used to select which groups of tests you want in occTest workflow
#' @details See occTest::showTests for further information on tests used in the packages
#' @param geo logical. Should this family of tests be performed?
#' @param lu logical. Should this family of tests be performed?
#' @param env logical. Should this family of tests be performed?
#' @param time logical. Should this family of tests be performed?
#' @return list
#' @keywords user
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
setTestBlocks      <- function (geo = T,
                                lu = T,
                                env = T,
                                time = T){
  
  #for testing
  #geo = T; lu = T; env = T; time = T
  
  objName = ls()
  paramsDF = list ()
  for (i in 1:length(objName)){
    paramsDF[[i]] = data.frame (testBlock = objName[i] ,activate = get(objName[i]))
  }
  
  paramsDF = do.call(rbind,paramsDF)
  allMetadata = readRDS(system.file('ext/fieldmetadata.rds',package='occTest'))
  newSettings = dplyr::left_join(allMetadata,paramsDF,by='testBlock')
  
  newParamsList = setTestTypes(countryStatusRange = unique (newSettings$activate [which (newSettings$testType == 'countryStatusRange')]),
                           centroidDetection = unique (newSettings$activate [which (newSettings$testType == 'centroidDetection')]),
                           humanDetection = unique (newSettings$activate [which (newSettings$testType == 'humanDetection')]),
                           landUseType = unique (newSettings$activate [which (newSettings$testType == 'landUseType')]),
                           institutionLocality = unique (newSettings$activate [which (newSettings$testType == 'institutionLocality')]),
                           geoOutliers = unique (newSettings$activate [which (newSettings$testType == 'geoOutliers')]),
                           envOutliers = unique (newSettings$activate [which (newSettings$testType == 'envOutliers')]),
                           geoenvLowAccuracy = unique (newSettings$activate [which (newSettings$testType == 'geoenvLowAccuracy')]))

  
    
  return (newParamsList)
}

# showTests ====
#' @title Show implemented tests and types of tests 
#' @description prints a table with the column names
#' @details The function prints a guide to column naming conventions used by occTest in their default parameters. These defaults can be changed via setTableNames, but the user may also decide to format their input table according to these naming conventions. 
#' @return prints a dataframe
#' @keywords user
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

showTests<- function (){
  tabNames=readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
  DT::datatable(tabNames)
}


# minimalSettings ====
#' @title Load minimal settings for occTest 
#' @description Loads a list of lists with the different default parameters for analysis, outputs and grading needed in occTest. It avoids using  some functions of the pkg under development.
#' @details it can be use internally or it can be used by a user to subsequently modify parameters
#' @return list of lists with all different parameters to use in occTest function
#' @keywords user
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

minimalSettings <- function (){
  
  #download info
  dest_url_hii = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/hii.rds'
  outFile_hii = paste0(tempdir(),'/hii.rds')
  utils::download.file(url=dest_url_hii,destfile = outFile_hii)
  
  defaultSettings = list (
    
    #grading Settings
    # gradingSettings = list (grading.test.type = 'majority', #other options are 'strict' 'relaxed'
    #                         qualifiers=T,
    #                         qualifier.label.scoping=c('A','B','C','D','E'))
    # 
    # ,
    #writing outputs settings
    writeoutSettings = list (#writing outputs
      output.dir=NULL,
      writeAllOutput=F, #overwrites write.simple.output, write.full.output
      write.simple.output=F,
      write.full.output=F,
      output.base.filename="occTest")
    ,
    #tableSettings
    tableSettings = list (taxonobservation.id = 'taxonobservationID',
                          x.field = 'decimalLongitude',
                          y.field = 'decimalLatitude',
                          t.field = 'eventDate', #time field (date)  
                          l.field = 'verbatimLocality', #locality field 
                          c.field = 'countryCode', #country field field (date)   
                          e.field = 'recordedElevation', #elevation recoreded  in meters
                          a.field = 'coordinateUncertaintyInMeters', #coordinate uncertainty in meters
                          ds.field = 'datasetName' #dataset field identifier
    )
    ,                     
    #analysis settings
    analysisSettings =list (
      geoSettings = list (
        coordinate.decimal.precision = 4,
        points.proj4string=sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      )
      ,
      countryStatusRange = list (
        countries.shapefile=rnaturalearth::ne_countries(scale = 50),
        countryfield.shapefile = 'iso_a3',
        doRangeAnalysis=T,
        excludeUnknownRanges= F,
        excludeNotmatchCountry= F,
        doCountryRecordAnalysis=T
      )
      ,
      centroidDetection = list (doCentroidDetection=T,
                                methodCentroidDetection='CoordinateCleaner'
      )
      ,
      humanDetection= list (doHumanDetection=T,
                            methodHumanDetection='all',
                            th.human.influence = 45,
                            ras.hii=readRDS(outFile_hii)
      )
      ,
      landUseType   = list (doLandUse=T,
                            methodLandUse='in',
                            landUseCodes = NULL,
                            ras.landUse=NULL #we need a default here to be downloaded
      )
      ,
      institutionLocality = list (doInstitutionLocality=T,
                                  methodInstitutionLocality='all'
      )
      ,
      
      geoOutliers = list (
        doGeoOutliers=T,
        methodGeoOutliers='all',
        alpha.parameter = 2
      )
      ,
      
      envOutliers = list (
        doEnvOutliers=T,
        methodEnvOutliers='all',
        th.perc.outenv =  0.2
      )
      ,
      geoenvLowAccuracy = list (methodGeoEnvAccuracy='all',
                                doGeoEnvAccuracy=T,
                                elev.quality.threshold = 100
      )
      
      
      
    )#end analysis settings
    
  )#end all default settings
  
  
  return (defaultSettings)
  
  
  
}



