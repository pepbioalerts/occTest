context("Workflow test")
require(terra)
require(sf)
require (occTest)

# setup test data
occData <- read.csv (system.file('ext/exampleOccData.csv',package = 'occTest'))
renv <- terra::rast (system.file('ext/AllEnv.tif',package = 'occTest'))
dem <- terra::rast (system.file('ext/DEM.tif',package = 'occTest'))
settings <- readRDS (system.file('ext/exSettings.rds',package = 'occTest'))



out = occTest(sp.name='MyFake species',
             sp.table = occData,ntv.ctry = 'ESP',inv.ctry = 'FRA',
             tableSettings = settings$tableSettings,
             writeoutSettings = settings$writeoutSettings,
             analysisSettings = settings$analysisSettings,
             r.env = renv,r.dem=dem)

tests_performed <- colnames(out)
out_docs <- readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
tests_documented <- out_docs %>% 
  tidyr::unite(documented_tests,testType:method) %>% 
  dplyr::pull (documented_tests) %>% paste0(c('_value','_test'))

tests_documented[!tests_documented %in% names(out)]

# run tests workflow inputs and outputs
test_that("occTest inputs vs outputs", {
  expect_equal(nrow (out), nrow (occData))
  expect_equal( c('occTest','data.frame'),class(out))
  expect_in(names(occData),names(out))
})

# expected output columns



