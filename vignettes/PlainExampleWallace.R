#load data
library(spocc)
df <- occ(query = 'Martes martes')
occ.data <-occ2df(df)

#we change the names here but we could just change it in the tableSettings
names(occ.data)[2:3] <- c('x','y')

#load needed environmental data
library(raster)
environmentRaster = raster::getData(name='worldclim',var='bio', res=10)

#start testing  occurrences
library(occTest)
outMartesMartes = occurrenceTests (sp.name='Martes martes', sp.table = occ.data,r.env = environmentRaster)

#run wrapper function (test + filter together)
outMartesMartes2 = occSimpFilter (spOcc = occ.data,env = environmentRaster)

#for wallace you may check occSimpFilter and implement it in two steps
# in occSimpFilter the names of the columns have defaults, and you can change them for x, y, data, isocountry.
# the date methods, however, have not yet been implemented fully nor the land use
# the parameter classification can also be a Wallace tick box where the user specifies a how strict they want to be
# as with everything this implementation is very very LIGTH, the more we populate the function with column names of things the more we can do about it. 


