###################################
# SIMPLE EXAMPLE
###################################

#load data
library (spocc)
df=readRDS(system.file('ext/pine_marten.RDS',package='occTest'))
occ.data <-spocc::occ2df(df)

#(down)load needed environmental data
environmentRaster = raster::getData(name='worldclim',var='bio', res=10,path = '~/RS/')#tempdir()

### STEP 1 TABLE FORMATING 
library(occTest)
# you can adapt the format names (adapt data frame): start formatting to the names in occTest
showTableNames ()

#in this example we only have three fields
occ.data.Newnames = occ.data
names (occ.data.Newnames)[2] <- 'decimalLongitude'
names (occ.data.Newnames)[3] <- 'decimalLatitude'
names (occ.data.Newnames)[4] <- 'eventDate'

#alternatively you can modify the settings to conform the names to your format
#for that, however you need to change the default parameters
#you can use showTableNames() to look at what do the parameters mean
names (occ.data)
showTableNames ()
myTableNames = setTableNames(x.field='longitude',y.field='latitude',t.field = 'date')

### STEP 2 SELECT ANALYSIS
#you can activate or deactivate certain analysis by modifying the analysis parameters
#see list of tests:
showTests ()

#using the setTests function you can deactivate certain types of tests
mySelectedAnalysis  = setTestTypes (centroidDetection = F,geoenvLowAccuracy = F)

#using the setTestBlocks function you can deactivate certain blocks (e.g. a kind of test types)
mySelectedAnalysis2 = setTestBlocks(time = F)

### STEP 3 RUN TESTS
martesTests = occTest(sp.name = "Martes_martes",sp.table = occ.data,r.env = environmentRaster,
                       #default parameters of column names changed
                       tableSettings = myTableNames , 
                       #default parameters of the analysis changed
                       analysisSettings = mySelectedAnalysis 
                       )


#STEP 4 SCRUB (FILTER) OCCURRENCES
occMartesFiltered = occTest::occFilter(df = martesTests,
                              by = 'testBlock',# no need this is default other option is testTypes
                              errorAcceptance = 'relaxed')

occMartesFiltered = occTest::occFilter(df = martesTests,
                                       by = 'testBlock',# no need this is default other option is testTypes
                                       errorAcceptance = 'strict')
#the output table filtered
occMartesFiltered$fitleredDataset

#the output scores for all records (except for those automatically filtered because of coordIssues)
occMartesFiltered$summaryStats



