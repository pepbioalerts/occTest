# occTest v0.1.2
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/ijtiff/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ijtiff/actions)


occTest is a workflow for general testing of species occurrence records in geographical, environmental and temporal scale. The workflow is structured through set of algorithms to identify potential problems with species occurrence records by employing a hierarchical organization of multiple tests. 

The workflow has a hierarchical structure organized in testPhases (i.e. cleaning vs. testing), that encompass different testBlocks grouping different testTypes (e.g. environmental outlier detection) which may use different testMethods (e.g. Rosner test, jacknife,etc.). Four different testBlocks characterize potential problems in geographic, environmental, human influence, and temporal dimensions.  Filtering and plotting functions are incorporated to facilitate the interpretation of tests. 

You can read the full details in this publication, which we encourage you to cite when using occTest:

Serra‐Diaz, Josep M., Jeremy Borderieux, Brian Maitner, Coline CF Boonman, Daniel Park, Wen‐Yong Guo, Arnaud Callebaut, Brian J. Enquist, Jens‐C. Svenning, and Cory Merow. "**occTest: An integrated approach for quality control of species occurrence data**" Global Ecology and Biogeography (2024): e13847. DOI: [https://doi.org/10.1111/geb.13847](url)


# Installation
## Stable from CRAN
Under review

## Development version from GitHub
```r
devtools::install_github("pepbioalerts/occTest")
library(occTest)
```

# Usage
Simple examples of a workflow below. A full explained vignette wit all options explained can be found at github.com/pepbioalerts/occtest. 

This second repo also serves as a data repo for automatic datasets used by occTest, such as human influence index or county centroids, when the user decides to use certain tests.

## Example 1: simple automatic testing
This is first automatic pass to occurrence testing with many default settings. These settings may or may not be appropriate to your analysis, so inspect the results carefully (e.g. plot them).

```r
#Download occurrence data
occ_data <- spocc::occ(query = 'Martes martes') 
occ_data <- spocc::occ2df(occ_data)
#Set the table columns to match occTest namings
names (occ_data)[c(2,3)] <- c('decimalLongitude','decimalLatitude')
#Load needed environmental data
environmentRaster <- geodata::worldclim_global(var = 'bio',path=tempdir(),res=10)
#Test occurrences (with default parameters)
outMartesMartes <- occTest::occTest(sp.name='Martes martes',
                                    sp.table = occ_data,
                                    r.env = environmentRaster)
class(outMartesMartes)
tail(names(outMartesMartes),n=7)
table(outMartesMartes$geoOutliers_score,exclude = F)
```
## Example 2: User-defined inputs
We use a predefined example dataset stored in the occTest package to showcase customization of the default parameters. 
```r
### Example 2 ====
#load data
sp_file <- system.file('ext/exampleOccData.csv',package='occTest')
occ_species_raw <- read.table(file = sp_file,header = T,sep = ',',as.is = T)
#load the raster of the environmental variables to the specific region
ras_env <- system.file('ext/AllEnv.tif',package='occTest') |> terra::rast()
#load a high resolution raster of elevation
ras_dem <- system.file('ext/DEM.tif',package='occTest') |> terra::rast()
#load a high resolution coastlines
pol_ctry <- system.file('ext/countries.rds',package='occTest') |> readRDS()
#show default settings
show_tableSettings()
#modify settings
my_table_Settings <- set_tableNames(x.field = 'MAPX',
                                    y.field = 'MAPY',
                                    t.field = 'DATE',
                                    l.field = 'LOCALITYNAME',
                                    c.field = 'CONTRYRECORD',
                                    e.field = 'ELEVATION',
                                    a.field =c('UNCERTAINTY_X_M','UNCERTAINTY_Y_M'))
my_analysis_settings <- set_testTypes(geoenvLowAccuracy = T)
my_writeout_settings <- set_writeout(output.dir = getwd(),writeAllOutput = T)
#run workflow
occ_species_test <- occTest(sp.name = 'Genus_species',
                            sp.table = occ_species_raw,
                            r.env =ras_env,
                            tableSettings = my_table_Settings,
                            analysisSettings = my_analysis_settings,
                            writeoutSettings = my_writeout_settings,
                            r.dem =   ras_dem,
                            ntv.ctry = 'ESP',
                            inv.ctry = 'FRA',
                            verbose = T)
#load default structure of settings
default_settings <- defaultSettings()
my_analysis_settings <- default_settings$analysisSettings
#modify some analysis table parameters
my_analysis_settings$humanAnalysis$th.human.influence <-50
my_analysis_settings$rangeAnalysis$countries.shapefile <-pol_ctry
my_analysis_settings$rangeAnalysis$countryfield.shapefile <- 'ISO'
my_analysis_settings$geoenvLowAccuracy$doGeoEnvAccuracy <- F
my_analysis_settings$envOutliers$methodEnvOutliers <- "grubbs"
#run occTest with my new parameters
occ_species_test <- occTest(sp.name = 'Genus_species',
                            sp.table = occ_species_raw,
                            r.env =ras_env,
                            tableSettings = my_table_Settings,
                            analysisSettings = my_analysis_settings,
                            writeoutSettings = my_writeout_settings,
                            r.dem =   ras_dem,
                            ntv.ctry = 'ESP',
                            inv.ctry = 'FRA',
                            verbose = T)
```
## Example 3: Filterig data
We provide here some example code to use simple rules for occurrence filtering based on consensus among methods. Two parameters control how occurrences are removed: (1) the error acceptance threshold (errorThreshold) and (2) the level at which the threshold is applied (by, e.g. testTypes or testBlocks). errorThreshold is a numeric value from 0-1 removing records above a certain specified level which can be at the testType level or the testBlock level. 
We strongly encourage to use these methods as a first check but it is recommended to investingate the test results and apply the rules that make more sense to your study.

```r
occ_species_test <- occTest(sp.name = 'Genus_species',
                            sp.table = occ_species_raw,
                            r.env =ras_env,
                            tableSettings = my_table_Settings,
                            r.dem =   ras_dem,
                            ntv.ctry = 'ESP',
                            inv.ctry = 'FRA',
                            verbose = T)
outFiltered_byTestType  <- occFilter(df = occ_species_test,errorThreshold = 0.8, 
                                     by = 'testType' )
str(outFiltered_byTestType,max.level = 1)
```
## Example 4: Plotting occTest results
We provide a generic plot function for occTest. It can easily be used to display the results of the occTest function. Additionally, if the results of the occFilter functions are provided, it will show on a map the results of the filtering rules.  
```r
#getting tree occurrences from a central European oak tree
df_qupe <- spocc::occ (query ='Quercus petraea',limit = 1000)
occ_data_qupe <-spocc::occ2df(df_qupe)
names (occ_data_qupe)[c(2,3)] <- c('decimalLongitude','decimalLatitude')
#downloading environmental data
renv <- geodata::worldclim_global(var = 'bio',path=tempdir(),res=10) 
#running occTest
out_qupe <- occTest(sp.name='Quercus petraea',
                    sp.table = occ_data_qupe ,
                    r.env = renv,
                    verbose = F)
#plot results
plot(out_qupe)
# running a strict filtering process
out_filtered_qupe <- occFilter(df = out_qupe , errorAcceptance = 'strict')
# adding the occFilter object to the plot
list_of_plots <- plot(x = out_qupe , occFilter_list = out_filtered_qupe, show_plot = F)
list_of_plots
```
## Example 5: Addding your own tests
You can add your own tests to an occTest functions, and reTest to get new ensembles of tests.

```r
myNewTest <- data.frame (humanDetection_myInventedTest_test=sample (c(T,F),replace = T,size = nrow(out_qupe)))
out_test_newAdded <- reTest(occTest_result = out_qupe, my_new_test =  myNewTest)
df_compare <- data.frame (new_score= out_test_newAdded$humanDetection_score ,old_score = out_qupe$humanDetection_score)
head (df_compare)
```

# Documentation
Further read at:
[https://github.com/pepbioalerts/vignetteXTRA-occTest](https://github.com/pepbioalerts/vignetteXTRA-occTest)

To download it (you can't just right-click on github), its easiest to install the chrome extension Gitzip:
[https://chrome.google.com/webstore/detail/gitzip-for-github/ffabmkklhbepgcgfonabamgnfafbdlkn](https://chrome.google.com/webstore/detail/gitzip-for-github/ffabmkklhbepgcgfonabamgnfafbdlkn)
You can then right-click on the vignette html to download 

# Citation
Read the full paper here:
Serra‐Diaz, J. M., Borderieux, J., Maitner, B., Boonman, C. C., Park, D., Guo, W. Y., ... & Merow, C. (2024). occTest: An integrated approach for quality control of species occurrence data. Global Ecology and Biogeography, e13847.
[https://doi.org/10.1111/geb.13847](https://doi.org/10.1111/geb.13847)
