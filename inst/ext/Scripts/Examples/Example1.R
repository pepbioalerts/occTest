#you should probably change this
setwd("D:/Test_functionWorkflow")


#### load needed data to run the function
require(raster)
occ.species.nogeoscurb <- read.csv ('./Example1/Pinus_halepensis_aggregated_PO_NOGeoScrub.csv')
environmental.variables.files <- paste0("./Environment/wc2.0_30s_bio/wc2.0_bio_30s_" , c("01","05","06","12","13","14"), ".tif" )
raster.environement <- raster::stack (as.list(environmental.variables.files))
raster.DEM <- raster ("./Altitude/Altitude.tif")

require(rgeos)
require(rgdal)
countries.pol <- readOGR(dsn ="D:/Test_functionWorkflow/CountryLevel",layer = 'Countries_SPDF' )
raster.humaninfluence <- raster ('./Environment/HII/hii_wgs84.tif')


####Run function
source ("./Scripts/F_workflow.R")
piha.occ.scrubbed <- occurrence.profiler(output.dir ="D:/Test_functionWorkflow/Example1/Output",
                                         sp.table = occ.species.nogeoscurb,
                                         sp.name = "Pinus halepensis",
                                         r.env = raster.environement,
                                         r.dem = raster.DEM ,
                                         countries.shapefile = countries.pol,
                                         countryfield.shapefile = 'ISO',
                                         ntv.ctry = c('ESP','ITA'),
                                         inv.ctry = NA ,
                                         ras.hii = raster.humaninfluence,
                                         x.field = 'x' ,
                                         y.field = 'y',
                                         t.field =  NULL,
                                         c.field = 'country_matched_iso3c',
                                         l.field = 'locality',
                                         e.field = 'elev',
                                         coordinate.decimal.precision = 2,
                                         points.proj4string = sp:::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                                         alpha.parameter = 2,
                                         th.human.influence = 45,
                                         th.perc.outenv = 0.2,
                                         elev.quality.threshold =250,
                                         qualifiers=T,
                                         write.simple.output = T,
                                         write.full.output = T,
                                         output.base.filename = "QAQC")


#### Inspect results
#full data.frame of the profiling with all the coments
head(piha.occ.scrubbed$occ_full_profile)

#short data.frame of the profiling with all the coments
head(piha.occ.scrubbed$occ_short_profile)

#take a look at the diferent tags generated
table(piha.occ.scrubbed$occ_short_profile$tag)

#take a look at the diferent grades generated
table(piha.occ.scrubbed$occ_short_profile$quality.grade)

#take a look at the diferent profiles generated
table(piha.occ.scrubbed$occ_short_profile$quality.label)




