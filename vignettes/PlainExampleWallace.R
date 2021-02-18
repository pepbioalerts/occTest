#load data
library(spocc)
#df <- occ(query = 'Martes martes')
df=readRDS(system.file('ext/pine_marten.RDS',package='occTest'))
occ.data <-occ2df(df)

#(down)load needed environmental data
environmentRaster = raster::getData(name='worldclim',var='bio', res=10,path = tempdir())

#start testing  occurrences
library(occTest)
tableSettings=occTest::defaultSettings()$tableSettings
tableSettings$x.field='longitude'
tableSettings$y.field='latitude'
outMartesMartes = occurrenceTests(sp.name='Martes martes', sp.table = occ.data,
                                  r.env = environmentRaster,
                                  tableSettings = tableSettings)
out2=occFilter(df = outMartesMartes,level = 1,errorAcceptance = 'majority')


#run wrapper function (test + filter together)
# CM: this one won't work with the demo unless you match the occ.data colnames to the standardized colnames of occTest
#outMartesMartes2 = occSimpFilter(spOcc = occ.data,env = environmentRaster)

#for wallace you may check occSimpFilter and implement it in two steps
# in occSimpFilter the names of the columns have defaults, and you can change them for x, y, data, isocountry.
# the date methods, however, have not yet been implemented fully nor the land use
# the parameter classification can also be a Wallace tick box where the user specifies a how strict they want to be
# as with everything this implementation is very very LIGTH, the more we populate the function with column names of things the more we can do about it. 



