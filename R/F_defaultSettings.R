#' @title load default settings for occTest
#'
#' @description Loads a list of lists with the different default parameters for analysis, outputs and grading needed in occTest
#' @details it can be use internally or it can be used by a user to subsequently modify parameters
#' @return list of lists with all different parameters to use in occProfile function
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


defaultSettings <- function (x){
  
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
    tableSettings = list (taxonobservation.id = 'taxonID',
                          x.field = 'decimalLatitude',
                          y.field = 'decimalLongitude',
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
        points.proj4string=sp:::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      )
      ,
      rangeAnalysis = list (
        countries.shapefile=rnaturalearth::ne_countries(),
        countryfield.shapefile = 'iso_a3',
        doRangeAnalysis=T,
        excludeUnknownRanges= F,
        excludeNotmatchCountry= F,
        doCountryRecordAnalysis=T
      )
      ,
      centroidAnalysis = list (doCentroidDetection=T,
                               methodCentroidDetection='all'
      )
      ,
      humanAnalysis= list (doHumanDetection=T,
                           methodHumanDetection='all',
                           th.human.influence = 45,
                           ras.hii=raster::raster(system.file('ext/hii/hii_wgs84_v2.tif',package='occTest'))
      )
      ,
      landUseAnalysis= list (doLandUse=T,
                             methodLandUse='in',
                             landUseCodes = NULL,
                             ras.landUse=NULL #we need a default here to be downloaded
      )
      ,
      institutionAnalysis = list (doInstitutionLocality=T,
                                  methodInstitutionLocality='all'
      )
      ,
      
      geooutliersAnalysis = list (
        doGeoOutliers=T,
        methodGeoOutliers='all',
        alpha.parameter = 2
      )
      ,
      
      envoutliersAnalysis = list (
        doEnvOutliers=T,
        methodEnvOutliers='all',
        th.perc.outenv =  0.2
      )
      ,
      accuracyAnalysis = list (methodGeoEnvAccuracy='all',
                               do.geoEnvAccuracy=T,
                               elev.quality.threshold = 100
      )
      
      
      
    )#end analysis settings
    
  )#end all default settings
  
  
  return (defaultSettings)
  
  
  
}
