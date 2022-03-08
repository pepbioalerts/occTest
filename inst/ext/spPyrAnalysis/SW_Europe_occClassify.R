#analysis
library (devtools)
library (BIEN)
library (raster)
library (data.table)
library (bit64)
library (rgbif)
library (parallel)
library (profvis)

#use the downloaded data from GBIF 
occGBIF = fread('/Users/pserra/RS/cleanOccAnalysis/GBIFdownolad/OccPyreneesArea.csv')

#get species 
spPyr = unique (occGBIF$species)
toEliminate = which(spPyr=='')
spPyr = spPyr [(-1) * toEliminate]

#Histogram of occurrences by species
Nsp =  (table (occGBIF$species))
hist (Nsp)
summary  (as.numeric (Nsp))
rm (occGBIF)
gc()
save(list = c('spPyr'), file = '/Users/pserra/RS/cleanOccAnalysis/GBIFdownolad/SpPyr.RData')

#get occurrence counts among all gbif data (occ counts with geographic data)

    # NspGBIF = lapply (spPyr, function (x){
    #   print(which (spPyr==x))
    #   key <- name_suggest(q=x, rank='species')$key[1]
    #   Ngeoref <- occ_count(taxonKey = key,georeferenced = T)
    #   data.frame (sp=x,Ngeoref=Ngeoref)  
    # })
    # 
    # NspGBIF = do.call (rbind,NspGBIF)
    # spidLowN = which (NspGBIF$Ngeoref<7)
    # length(spidLowN)
    # spidMedN = which (NspGBIF$Ngeoref<200 & NspGBIF$Ngeoref>=7)
    # length(spidMedN)
    # spidHighN = which (NspGBIF$Ngeoref<2000 & NspGBIF$Ngeoref>=200)
    # length(spidHighN)
    # spidSuperHighN = which (NspGBIF$Ngeoref>=2000)
    # length(spidHighN)

#get climate 
lf = list.files ("/Users/pserra/RS/GIS/CHELSA/",full.names = T)
env = raster::stack (lf)

#get digital elevation model
elev = raster::raster('/Users/pserra/RS/GIS/mn30_grd/mn30.tif')

#mySettings
mySettings <- occTest::defaultSettings()
mySettings$tableSettings$x.field <- 'longitude'
mySettings$tableSettings$y.field <- 'latitude'
mySettings$tableSettings$t.field <- 'date_collected'
mySettings$tableSettings$l.field <- 'locality'
mySettings$tableSettings$e.field <- 'elevation_m'
#I DO NOT SEE IT IN BIEN
mySettings$tableSettings$a.field <- c('UNCERTAINTY_X_M','UNCERTAINTY_Y_M')


#download and save the data for each species 
spPyrDone = list.dirs(path ="/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/" ,full.names = F,recursive = F)
spPyrNotDone = setdiff(spPyr,spPyrDone)

out = lapply (spPyrNotDone, function (spModel){
  print (spModel)
    outDownload  = try({
    
    print(which (spPyr == spModel) )
    occRGBIF = occ_data(scientificName = spModel,limit = 10000)
    #need to go from 2-digit to three digit for occ profiler
    occRGBIF$data$countryCode3c = countrycode::countrycode (occRGBIF$data$countryCode,origin ='iso2c',destination = 'iso3c')
    #get the citations for the species
    cit = gbif_citation(occRGBIF)
    #write inputs 
    pathInputs = paste0('/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/',spModel)
    jspModel = occTest:::.join.spname(spModel)
    dir.create (pathInputs,showWarnings = F,recursive = T)
    save(occRGBIF,file = paste0(pathInputs,'/', jspModel,'.RData'))
    save(cit,file = paste0(pathInputs,'/', jspModel,'_citations.RData'))  },silent = T)
    
    
    if (class(outDownload)=="try-error") {return ("Error")} else {return ('Download ended')}
    
    
    

})



#this species could not be downloaded by rgbif, we use BIEN
spPyrNotDoneRGBIF = spPyrNotDone
lapply (13:length(spPyrNotDoneRGBIF), function (i){
  print (i)
  spModel= spPyrNotDoneRGBIF[i]
  occBIEN = BIEN_occurrence_species(species = spModel,cultivated = T,only.new.world = F,native.status = T,political.boundaries = T,collection.info = T)
  if (nrow (occBIEN)==0) {return (NULL)}
  #need to go from 2-digit to three digit for occ profiler
  occBIEN$country = countrycode::countrycode (occBIEN$country ,origin ='country.name',destination = 'iso3c')
  #write inputs 
  jspModel = occTest:::.join.spname(spModel)
  pathInputs = paste0('/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/',jspModel)
  dir.create (pathInputs,showWarnings = F,recursive = T)
  save(occBIEN,file = paste0(pathInputs,'/', jspModel,'_BIEN.RData'))
  
  
})

#recheck folder names
pathInputs = "/Users/pserra/RS/cleanOccAnalysis/QAQCinputs"
ldirs = list.dirs(path = pathInputs,full.names = T,recursive = F)
lapply (ldirs, function (x){
  ssplit = strsplit (x, ' ')[[1]]
  if (length(ssplit)) {
    ssplpit2 = paste(ssplit,collapse = '_')
    file.rename(from = x,to = ssplpit2)
  }
})


#compute QAQC in parallel for those species with >500 records
out = lapply (spPyr[8:length(spPyr)], function (spModel){
  try({
    
    print(which (spPyr == spModel) )
    
    occRGBIF = occ_data(scientificName = spModel,limit = 10000)
    #need to go from 2-digit to three digit for occ profiler
    occRGBIF$data$countryCode3c = countrycode::countrycode (occRGBIF$data$countryCode,origin ='iso2c',destination = 'iso3c')
    #get the citations for the species
    cit = gbif_citation(occRGBIF)
    #write inputs 
    pathInputs = paste0('/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/',spModel)
    jspModel = occTest:::.join.spname(spModel)
    dir.create (pathInputs,showWarnings = F,recursive = T)
    save(occRGBIF,file = paste0(pathInputs,'/', jspModel,'.RData'))
    save(cit,file = paste0(pathInputs,'/', jspModel,'_citations.RData'))
    
    if (nrow (occRGBIF$data) <250) {return(spModel)}
    
    
    #build the Settings for QAQC
    rgbifSettings <- occTest:::defaultSettings()
    rgbifSettings$tableSettings$x.field <- 'decimalLongitude'
    rgbifSettings$tableSettings$y.field <- 'decimalLatitude'
    rgbifSettings$tableSettings$t.field <- 'eventDate'
    rgbifSettings$tableSettings$l.field <- 'locality'
    rgbifSettings$tableSettings$e.field <- 'elevation' #or should we use verbatimElevation?
    rgbifSettings$tableSettings$c.field <- 'countryCode3c'
    rgbifSettings$tableSettings$a.field <- c('coordinateUncertaintyInMeters','coordinateUncertaintyInMeters')
    
    myOutputs = '/Users/pserra/RS/cleanOccAnalysis/QAQCoutputs'
    dir.create (myOutputs, recursive = T,showWarnings = F)
    
    rgbifSettings$writeoutSettings$output.dir = myOutputs
    rgbifSettings$writeoutSettings$write.simple.output = T
    rgbifSettings$writeoutSettings$write.full.output = T
    
    spProfile_1 = occurrenceClassify(sp.name= spModel,
                                     sp.table = occRGBIF$data, 
                                     r.env = env , 
                                     r.dem = elev,
                                     tableSettings =rgbifSettings$tableSettings,
                                     writeoutSettings = rgbifSettings$writeoutSettings,
                                     doParallel = T, mc.cores = 6,
                                     verbose = T )
    

    if (class(spProfile_1)=="try-error") {return ("Error")} else {return ('Model ended')}
    
    
    
  },silent = T)
})


#### PROFILING THE FUNCTION 
# p2 <- profvis ({spProfile_1 = try (occurrenceClassify(sp.name= spModel,
#                                                      sp.table = occRGBIF$data, 
#                                                      r.env = env , 
#                                                      r.dem = elev,
#                                                      tableSettings =rgbifSettings$tableSettings,
#                                                      writeoutSettings = rgbifSettings$writeoutSettings), silent=T)})
# 
# 

# we deactivate grade E
#rgbifSettings$gradingSettings$qualifier.label.scoping = c("A","B","C" ,"D")
# spProfile_3 = occurrenceClassify(sp.name= spModel,
#                                  sp.table = occRGBIF$data, 
#                                  r.env = env , 
#                                  r.dem = elev,
#                                  tableSettings =rgbifSettings$tableSettings,
#                                  gradingSettings = rgbifSettings$gradingSettings,
#                                  interactiveMode = F,
#                                  resolveAlienCtry = F,resolveNativeCtry = F)
# 
# table(spProfile_3[[1]]$quality.label)
