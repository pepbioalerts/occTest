
require (rnaturalearth)


defaultSettings = list (
  
  

#grading Settings
gradingSettings = list (grading.test.type = 'majority', #other options are 'strict' 'relaxed'
                        qualifiers=T,
                        qualifier.label.scoping=c('A','B','C','D','E'))

,
#writing outputs settings
writeoutSettings = list (#writing outputs
  output.dir=NULL,
  writeAllOutput=F, #overwrites write.simple.output, write.full.output
  write.simple.output=F,
  write.full.output=F,
  output.base.filename="QAQC")
,
#tableSettings
tableSettings = list (taxonobservation.id = NULL,
                      x.field = 'x',
                      y.field = 'y',
                      t.field = NULL,
                      l.field = NULL,
                      c.field = NULL,
                      e.field = NULL,
                      a.field = NULL, #coordinate uncertainty in meters
                      ds.field = NULL #dataset field identifier
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
  humanAnalysis= list (doHyperHumanDetection=T,
                       methodHyperHumanDetection='all',
                       th.human.influence = 45,
                       ras.hii=raster::raster(system.file('ext/hii/hii_wgs84_v2.tif',package='occProfileR'))
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

