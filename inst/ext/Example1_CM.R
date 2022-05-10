#you should probably change this to your local directory
myDir='/Users/ctg/Dropbox/Projects/R_Packages/occTest/occTestData'

#### load needed data to run the function
#occ.species.nogeoscurb=read.csv(paste0(myDir,'/Example1/Pinus_halepensis_aggregated_PO_NOGeoScrub.csv'))
occ.species.nogeoscurb=spocc::occ2df(spocc::occ('Pinus_halepensis',from='gbif'))
environmental.variables.files <- paste0(myDir,"/Environment/wc2.0_30s_bio/wc2.0_bio_30s_" , c("01","05","06","12","13","14"), ".tif" )
raster.environement <- raster::stack (environmental.variables.files)
raster.DEM <- raster::raster(paste0(myDir,'/Altitude/Altitude.tif'))

countries.pol <- rgdal::readOGR(paste0(myDir,"/CountryLevel/Countries_SPDF.shp"))
raster.humaninfluence <- raster::raster(paste0(myDir,
                                               '/Environment/HII/hii_wgs84.tif'))


####Run function
piha.occ.scrubbed <- occurrence.profiler(output.dir=getwd(),
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
                                         points.proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
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




