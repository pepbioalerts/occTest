library (raster)
occCsv = read.csv ('~/RS/cleanOccAnalysis/issuesOccTest/Para_Pep_coordinate_reassignment.csv')

allEnv = raster::stack (list.files('~/RS/GIS/CLIMATE/CHELSA/',pattern = '.tif$',full.names = T))
rm (allEnv)
allEnv2 = terra::rast (list.files('~/RS/GIS/CLIMATE/CHELSA/',pattern = '.tif$',full.names = T))
kk = sum (allEnv2)


?terra
rasHabitat = rasallEnv
bio1 = raster::raster ('~/RS/GIS/CLIMATE/CHELSA/CHELSA_bio10_01.tif')
bio12 = raster::raster ('~/RS/GIS/CLIMATE/CHELSA/CHELSA_bio10_12.tif')
env = raster::stack (bio1,bio12)
mySettings = occTest::minimalSettings()
names (occCsv)
mySettings$tableSettings$x.field ='Longitude'
mySettings$tableSettings$y.field ='Latitude'
magu = occTest::occTest(sp.name='Matayba guianensis', sp.table = occCsv,r.env = env,
                                tableSettings = mySettings$tableSettings,
                                analysisSettings = mySettings$analysisSettings)


#investigate outputs via aggreagated scores per issue type
library (tidyverse)
magu %>% select (ends_with('_score')) %>% summary 
#investigate outputs via each test performed
magu %>% select (ends_with('_test')) %>% summary 


#check the consequences of applying different filters
relaxFilter = occTest::occFilter(magu,errorAcceptance = 'relaxed')
nrow (relaxFilter$fitleredDataset)
majFilter = occTest::occFilter(magu,errorAcceptance = 'majority') 
nrow (majFilter$fitleredDataset)
strictFilter = occTest::occFilter(magu,errorAcceptance = 'strict') 
nrow (strictFilter$fitleredDataset)

#plot the differences 
library (sf)
spMagu = sp::SpatialPointsDataFrame(coords = magu[,c('Longitude','Latitude')],data = magu)
#spMagu2 = sf::st_as_sf(magu,coords = c('Longitude','Latitude'))
bio1Crop = raster::crop (bio1,extent(spMagu2))
plot (bio1Crop)
plot (spMagu,add=T)
points (relaxFilter$fitleredDataset[,c('Longitude','Latitude')],col='red')
points (majFilter$fitleredDataset[,c('Longitude','Latitude')],col='purple')
points (strictFilter$fitleredDataset[,c('Longitude','Latitude')],col='black')

