#analysis
library (devtools)
library (BIEN)
library (raster)
library (data.table)
library (bit64)
library (rgbif)
library (parallel)
library (profvis)
library (parallelsugar)
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_

#use the downloaded data from GBIF 
#occGBIF = fread('/Users/pserra/RS/cleanOccAnalysis/GBIFdownolad/OccPyreneesArea.csv')

#get species 
spPyrTable = fread ('D:\\Write\\Pep\\PyrSDM\\SpeciesList\\GBIF_speciesList/0068556-200221144449610/0068556-200221144449610.csv')
str (spPyrTable)
spPyr = unique (spPyrTable$species)
toEliminate = which(spPyr=='')
spPyr = spPyr [(-1) * toEliminate]

#get gbif taxon key
spPyrTaxonKey = unique (spPyrTable$taxonKey)


# MASSIVE DOWNLOAD OF OCCURRENCES BY USERNAME AND PASSOWRD
# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!
occ_download(
  pred_in("taxonKey", spPyrTaxonKey),
  format = "SIMPLE_CSV",
  user='pep.agroparistech',pwd='Treechange2019@',email='pep.agroparistech@gmail.com'
  )


#Histogram of occurrences by species
  # Nsp =  (table (occGBIF$species))
  # hist (Nsp)
  # summary  (as.numeric (Nsp))
  # rm (occGBIF)
  # gc()
  # save(list = c('spPyr'), file = '/Users/pserra/RS/cleanOccAnalysis/GBIFdownolad/SpPyr.RData')

#get climate 

lf = paste0 ("D:/quercus/BD_SIG/climat/monde/Chelsa_1979_2013/climatologies/bio/CHELSA_bio10_",c(1,5,6,12,13,14),'.tif')
env = raster::stack (x = lf)


#get digital elevation model
elev = raster::raster('D:/quercus/BD_SIG/climat/monde/Chelsa_1979_2013/dem/mn30.tif')

#mySettings
mySettings <- occProfileR::defaultSettings()
mySettings$tableSettings$x.field <- 'longitude'
mySettings$tableSettings$y.field <- 'latitude'
mySettings$tableSettings$t.field <- 'date_collected'
mySettings$tableSettings$l.field <- 'locality'
mySettings$tableSettings$e.field <- 'elevation_m'

#I DO NOT SEE IT IN BIEN
mySettings$tableSettings$a.field <- c('UNCERTAINTY_X_M','UNCERTAINTY_Y_M')


#download and save the data for each species
library (foreach)
library (parallel)
library (doSNOW)

#create output directory and check which species have been done 
dir.create("D:/Write/Pep/PyrSDM/OccData")
spPyrDone = list.dirs(path ="D:/Write/Pep/PyrSDM/OccData" ,full.names = F,recursive = F)
source("D:\\Write\\Pep\\Github\\myHandyFuncs\\MyHandyFuncs.R")
spPyrHyphen = lapply (spPyr, join.spname) 
spPyrHyphen = unlist (spPyrHyphen)
spPyrNotDone = setdiff(spPyrHyphen,spPyrDone)


#============ old when we downloaded species one by one
### DOWNLOADS USING RGBIF AND RBIEN if NEEDBE


    # cl = makeCluster(5)
    # registerDoSNOW(cl)
    # iterations = length(spPyrNotDone)
    # pb = txtProgressBar(max = iterations, style = 3)
    # progress = function(n) setTxtProgressBar(pb, n)
    # opts = list(progress = progress)
    # lp = loaded_packages()
    # 
    # out = foreach::foreach(spModel = spPyrNotDone,.options.snow=opts, .errorhandling = 'remove',.packages = lp$package) %dopar% {
    #   
    #   outDownload  = try({
    #     print(which (spPyr == spModel) )
    #     occRGBIF = rgbif::occ_search(scientificName = spModel,limit = 10000)
    #     #get the citations for the species
    #     cit = gbif_citation(occRGBIF)
    #     #need to go from 2-digit to three digit for occ profiler
    #     occRGBIF$data$countryCode3c = countrycode::countrycode (occRGBIF$data$countryCode,origin ='iso2c',destination = 'iso3c')
    #     #write inputs 
    #     jspModel = occProfileR:::.join.spname(spModel)
    #     pathInputs = paste0('D:/Write/Pep/PyrSDM/OccData/',jspModel)
    #     dir.create (pathInputs,showWarnings = T,recursive = T)
    #     save(occRGBIF,file = paste0(pathInputs,'/', jspModel,'_gbif.RData'))
    #     save(cit,file = paste0(pathInputs,'/', jspModel,'_citations.RData'))  }
    #     ,silent = T)
    #   
    #   if (class(outDownload)=="try-error") {return ("Error")} else {return ('Download ended')}
    #   
    # }
    # registerDoSEQ()
    # stopCluster(cl)
    # spPyrDone = list.dirs(path ="D:/Write/Pep/PyrSDM/OccData" ,full.names = F,recursive = F)
    # jspPyr  =  unlist (lapply (spPyr, function (x) occProfileR:::.join.spname(x)))
    # spPyrNotDone = setdiff(jspPyr,spPyrDone)
    # idsSpNotDone = which (! jspPyr %in% spPyrDone)
    # spPyrNotDone = spPyr[idsSpNotDone]
    # 
    # 
    # #this species could not be downloaded by rgbif, we use BIEN
    # spPyrNotDoneRGBIF = spPyrNotDone
    # lapply (1:length(spPyrNotDoneRGBIF), function (i){
    #   print (i)
    #   spModel= spPyrNotDoneRGBIF[i]
    #   occBIEN = BIEN_occurrence_species(species = spModel,cultivated = T,only.new.world = F,native.status = T,political.boundaries = T,collection.info = T)
    #   if (nrow (occBIEN)==0) {return (NULL) ; print ('No records found')}
    #   #need to go from 2-digit to three digit for occ profiler
    #   occBIEN$country = countrycode::countrycode (occBIEN$country ,origin ='country.name',destination = 'iso3c')
    #   #write inputs 
    #   jspModel = occProfileR:::.join.spname(spModel)
    #   pathInputs = paste0('D:/Write/Pep/PyrSDM/OccData/',jspModel)
    #   dir.create (pathInputs,showWarnings = F,recursive = T)
    #   save(occBIEN,file = paste0(pathInputs,'/', jspModel,'_BIEN.RData'))
    #   
    #   
    # })
    # 
    # 
    # #check occurrence counts 
    # ldirs = list.dirs(path = pathInputs,full.names = T,recursive = F)
    # library (parallelsugar)
    # NoccSp= mclapply (ldirs,mc.cores= 10,FUN =  function (l){
    #   print (which (ldirs==l))
    #   
    #   allfiles = list.files(path = l,full.names = T)
    #   dfile = allfiles[1]
    #   spName = strsplit (basename(dfile),split = '.RData')[[1]]
    #   occDat = load(dfile)
    #   occDat = (get(occDat))
    #   if (class(occDat)=='gbif') {occDat = occDat$data; occinfo='RGBIF'} else {occinfo='BIEN'}
    #   Nocc=nrow(occDat)
    #   data.frame(spName, Nocc)
    #   
    # })
    # NoccSp = do.call (rbind, NoccSp)
    # 
    # 
    # #how many occurrence per species ?
    # Nclassification = cut (NoccSp$Nocc,c(10,50,200,500,100001))
    # table (Nclassification)
    # spLarge = as.character (NoccSp[which(NoccSp$Nocc>200),'spName'])
    # spAll = as.character (NoccSp$spName)

#============ 

#if massive download done then read in data

allSppfread ('D:/Write/Pep/PyrSDM/OccData/GBIForg_download/0068606-200221144449610/0068606-200221144449610.csv')

#################
#compute QAQC in occProfileR parallel MODE for those species with >200 records
#################

library(foreach)
library (parallel)
library (doSNOW)

# check species not done
pathOutputs = "D:/Write/Pep/PyrSDM/OccProfile"
ldirs = list.dirs(path = pathInputs, full.names = F,recursive = F)
ldirs = ldirs [grep(x = ldirs,'spatialData')* (-1)]
spNotDone =  spAll [which (! spAll %in% ldirs)]

#check cluster 
cl = makeCluster(5)
registerDoSNOW(cl)
iterations = length(spNotDone)
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
lp = loaded_packages()

out = foreach::foreach(spModel = spNotDone,.options.snow=opts, .errorhandling = 'remove',.packages = lp$package,.combine = c) %dopar% {
    
    #LOAD THE DATA 
    pathInputs = "D:/Write/Pep/PyrSDM/OccData"
    allfiles = list.files(path = paste0(pathInputs,'/',spModel) ,full.names = T)
    dfile = allfiles[1]
    spName = strsplit (basename(dfile),split = '.RData')[[1]]
    occDat = load(dfile)
    occDat = (get(occDat))
    if (class(occDat)=='gbif') {occDat = occDat$data; occinfo='RGBIF'} else {occinfo='BIEN'}
    
    
    #build the Settings for QAQC
    rgbifSettings <- occProfileR:::defaultSettings()
    rgbifSettings$tableSettings$x.field <- 'decimalLongitude'
    rgbifSettings$tableSettings$y.field <- 'decimalLatitude'
    rgbifSettings$tableSettings$t.field <- 'eventDate'
    rgbifSettings$tableSettings$l.field <- 'locality'
    rgbifSettings$tableSettings$e.field <- 'elevation' #or should we use verbatimElevation?
    rgbifSettings$tableSettings$c.field <- 'countryCode3c'
    if (occinfo=='RGBIF') {rgbifSettings$tableSettings$a.field <- c('coordinateUncertaintyInMeters','coordinateUncertaintyInMeters')}
    
    #PREPARE OUTPUTS
    myOutputs = "D:/Write/Pep/PyrSDM/OccProfile"
    dir.create (myOutputs, recursive = T,showWarnings = F)
    
    rgbifSettings$writeoutSettings$output.dir = myOutputs
    rgbifSettings$writeoutSettings$write.simple.output = T
    rgbifSettings$writeoutSettings$write.full.output = T
    
    #RUN FUNCTION
    spSci = occProfileR:::splitSpname(spModel)
    spProfile_1 = try (occProfileR:::occurrenceClassify(sp.name= spSci,
                                                   sp.table = occRGBIF$data, 
                                                   r.env = env , 
                                                   r.dem = elev,
                                                   tableSettings =rgbifSettings$tableSettings,
                                                   writeoutSettings = rgbifSettings$writeoutSettings,
                                                   doParallel = T, mc.cores = 6,
                                                   verbose = T ), silent=T)
    
    
    
    
  
  if (class(spProfile_1)=="try-error") {return ("Error")} else {return ('Success')}
  
  
}

registerDoSEQ()
stopCluster(cl)


#check species with errors
idSpError = which (out == 'Error')
spError = spNotDone [idSpError]

cl = makeCluster(6)
registerDoSNOW(cl)
iterations = length(spError)
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
lp = loaded_packages()

outErrors = foreach::foreach(spModel = spError,.options.snow=opts, .errorhandling = 'remove',.packages = lp$package,.combine = c) %dopar% {
  
  #LOAD THE DATA 
  pathInputs = "D:/Write/Pep/PyrSDM/OccData"
  allfiles = list.files(path = paste0(pathInputs,'/',spModel) ,full.names = T)
  dfile = allfiles[1]
  spName = strsplit (basename(dfile),split = '.RData')[[1]]
  occDat = load(dfile)
  occDat = (get(occDat))
  if (class(occDat)=='gbif') {occDat = occDat$data; occinfo='RGBIF'} else {occinfo='BIEN'}
  
  
  #build the Settings for QAQC
  rgbifSettings <- occProfileR:::defaultSettings()
  rgbifSettings$tableSettings$x.field <- 'decimalLongitude'
  rgbifSettings$tableSettings$y.field <- 'decimalLatitude'
  rgbifSettings$tableSettings$t.field <- 'eventDate'
  rgbifSettings$tableSettings$l.field <- 'locality'
  rgbifSettings$tableSettings$e.field <- 'elevation' #or should we use verbatimElevation?
  rgbifSettings$tableSettings$c.field <- 'countryCode3c'
  
  if (occinfo=='RGBIF') {rgbifSettings$tableSettings$a.field <- c('coordinateUncertaintyInMeters','coordinateUncertaintyInMeters')}
  
  #PREPARE OUTPUTS
  myOutputs = "D:/Write/Pep/PyrSDM/OccProfile"
  dir.create (myOutputs, recursive = T,showWarnings = F)
  
  rgbifSettings$writeoutSettings$output.dir = myOutputs
  rgbifSettings$writeoutSettings$write.simple.output = T
  rgbifSettings$writeoutSettings$write.full.output = T
  
  #RUN FUNCTION
  spSci = occProfileR:::splitSpname(spModel)
  spProfile_1 = try (occProfileR:::occurrenceClassify(sp.name= spSci,
                                                      sp.table = occRGBIF$data, 
                                                      r.env = env , 
                                                      r.dem = elev,
                                                      tableSettings =rgbifSettings$tableSettings,
                                                      writeoutSettings = rgbifSettings$writeoutSettings,
                                                      doParallel = T, mc.cores = 6,
                                                      verbose = T ), silent=T)
  
  
  
  
  
  if (class(spProfile_1)=="try-error") {return ("Error")} else {return (spProfile_1)}
  
  
}


registerDoSEQ()
stopCluster(cl)

##investigate which species have errors

outErrors[[5]]

