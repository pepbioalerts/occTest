

## UPDATE METHODS IN OCCTEST

#STEP1 INSERT METHOD IN THE APPRPRIATE TEST BLOCK

#UPDATE METADATA
library (tidyverse)

tabNames=readRDS(system.file('ext/analysisSettings_metadata.rds',package='occTest'))
View (tabNames)
names (tabNames)
tabNames2 <- tabNames %>% add_row(analysisSetting_level1='geoOutliers', 
                                  analysisSetting_level2='mcp_percSample', 
                                  definition='numeric. Percentage of samples used to build the minimum convex polygon', 
                                 .before=26)


saveRDS (tabNames2,'inst/ext/analysisSettings_metadata.rds')


fieldMetadata=readRDS('inst/ext/fieldMetadata.rds')
View (fieldMetadata)
names (fieldMetadata)
test
geo
geoOutliers
quantileSamplingCorr

fieldMetadata2 <- fieldMetadata %>% add_row(phase='test', 
                                  testBlock='geo', 
                                  testType='geoOutliers', 
                                  method='mcp',
                                  .before=22)

saveRDS (fieldMetadata2,'inst/ext/fieldMetadata.rds')
