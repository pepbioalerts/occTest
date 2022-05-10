#Example script that doesn't rely on data stored in folders

library(raster)
library(BIEN)
library(rgeos)
library(rgdal)
library(maptools)
library(biogeo)
library(speciesgeocodeR)
library(alphahull)
library(qdap)
library(plyr)

occ.species.nogeoscurb <- BIEN_occurrence_species(species = "Pinus contorta")
raster.environement <- getData('worldclim', var='bio', res=10)
raster.DEM <- getData('worldclim', var='alt', res=10)
data(wrld_simpl)
countries.pol <- wrld_simpl
rm(wrld_simpl)

#Why does code need iso code?

#Need to replace this with a link or something that avoids requiring a file
raster.humaninfluence <- raster ("C:/Users/Brian/Desktop/test_junk/occTest/Environment/HII/hii_wgs84_MED.tif")


source("GeoEnvironmentalProfileR/R/F_Analysis.R")
source("GeoEnvironmentalProfileR/R/F_Checks.R")
source("GeoEnvironmentalProfileR/R/F_Grading.R")
source("GeoEnvironmentalProfileR/R/F_OtherFuncs.R")
source("GeoEnvironmentalProfileR/R/F_qualifiers.R")
source("GeoEnvironmentalProfileR/R/F_workflow.R")
source("GeoEnvironmentalProfileR/R/Script_grading.R")


raster.humaninfluence@crs<-countries.pol@proj4string

raster.humaninfluence@crs
countries.pol@proj4string
raster.DEM@crs<-countries.pol@proj4string
raster.environement@crs<-countries.pol@proj4string

piha.occ.scrubbed <- occurrence.profiler(output.dir ="C:/Users/Brian/Desktop/test_junk/",
                                         sp.table = occ.species.nogeoscurb,
                                         sp.name = unique(occ.species.nogeoscurb$scrubbed_species_binomial),
                                         r.env = raster.environement,
                                         r.dem = raster.DEM ,
                                         countries.shapefile = countries.pol,
                                         countryfield.shapefile = 'ISO3',
                                         ntv.ctry = NULL,
                                         inv.ctry = NULL ,
                                         ras.hii = raster.humaninfluence,
                                         x.field = 'longitude' ,
                                         y.field = 'latitude',
                                         t.field =  NULL,
                                         c.field = NULL,
                                         l.field = NULL,
                                         e.field = NULL,
                                         coordinate.decimal.precision = 2,
                                         points.proj4string = countries.pol@proj4string,
                                         alpha.parameter = 2,
                                         th.human.influence = 45,
                                         th.perc.outenv = 0.2,
                                         elev.quality.threshold =250,
                                         qualifiers=T,
                                         write.simple.output = T,
                                         write.full.output = T,
                                         output.base.filename = "QAQC",
                                         qualifier.label.scoping = c("A"))

#cURRENT ISSUE:

#Error in eval(parse(text = grading.params$Conditional)) : 
#object 'missingEnvironment' not found 

#Seems to be associated with qgrade.data.frame() function in F_Grading.R file
