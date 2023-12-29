# occTest v0.1.2
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/ijtiff/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ijtiff/actions)

**occTest is currently under review and has been recently adapted to sf and terra packages. It will soon be submitted to CRAN once final checks are realized**

A workflow for general testing of species occurrence records in geographical, environmental and temporal scale. The workflow is structured through set of algorithms to identify potential problems with species occurrence records by employing a hierarchical organization of multiple tests. 

The workflow has a hierarchical structure organized in testPhases (i.e. cleaning vs. testing), that encompass different testBlocks grouping different testTypes (e.g. environmental outlier detection) which may use different testMethods (e.g. Rosner test, jacknife,etc.). Four different testBlocks characterize potential problems in geographic, environmental, human influence, and temporal dimensions.  Filtering and plotting functions are incorporated to facilitate the interpretation of tests. 


# Installation
## Stable from CRAN
Under review

## Developmental from GitHub
```r
devtools::install_github("pepbioalerts/occTest")
library(occTest)
```

# Usage
Simple example of a workflow below. A full explained vignette can be found at github.com/pepbioalerts/occtest. This second repo also serves as a data repo for automatic datasets used but occTest, such as human influence index or county centroids, when the user decides to use certain tests.

```r
### THIS IS A CUT DOWN  EXAMPLE 
### visit vignetteXtra-occTest for more info
#load environmental SpatRast

library (occTest)
#load occurrence data
occData <- read.csv (system.file('ext/exampleOccData.csv',package = 'occTest'))
#load environmental SpatRast
renv <- terra::rast (system.file('ext/AllEnv.tif',package = 'occTest'))
#load elevation SpatRast
dem <- terra::rast (system.file('ext/DEM.tif',package = 'occTest'))
#load settings
settings <- readRDS (system.file('ext/exSettings.rds',package = 'occTest'))
#run occTest
out = occTest(sp.name='MyFake species',
              sp.table = occData,ntv.ctry = 'ESP',inv.ctry = 'FRA',
              tableSettings = settings$tableSettings,
              writeoutSettings = settings$writeoutSettings,
              analysisSettings = settings$analysisSettings,
              r.env = renv,r.dem=dem)
```

# Documentation
[https://github.com/pepbioalerts/vignetteXTRA-occTest](https://github.com/pepbioalerts/vignetteXTRA-occTest)

To download it (you can't just right-click on github), its easiest to install the chrome extension Gitzip:
[https://chrome.google.com/webstore/detail/gitzip-for-github/ffabmkklhbepgcgfonabamgnfafbdlkn](https://chrome.google.com/webstore/detail/gitzip-for-github/ffabmkklhbepgcgfonabamgnfafbdlkn)
You can then right-click on the vignette html to download 

# Citation
Under review 
