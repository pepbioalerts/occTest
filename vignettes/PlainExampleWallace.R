#load data
library(spocc)
df <- occ(query = 'Martes martes')
occ.data <-occ2df(df)

#we change the names here but we could just change it in the tableSettings
names (occ.data)[2] <- 'x'
names (occ.data)[3] <- 'y'

#load needed environmental data
library (raster)
environmentRaster = raster::getData(name='worldclim',var='bio', res=10)

#start testing  occurrences
library(occProfileR)
outMartesMartes = occurrenceTests (sp.name='Martes martes', sp.table = occ.data,r.env = environmentRaster)
str (outMartesMartes)s
str(outMartesMartes$occTest_full)
#run wrapper function (test + classify, ont implemented yet, now it gives the same results)
outMartesMartes2 = occSimpleClassification (spOcc = occ.data,env = environmentRaster)


