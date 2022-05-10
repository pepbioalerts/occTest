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

user = Sys.info()['user']


#use the downloaded data from GBIF 
#occGBIF = fread('/Users/pserra/RS/cleanOccAnalysis/GBIFdownolad/OccPyreneesArea.csv')

#get species 
spPyr = fread ('D:\\Write\\Pep\\PyrSDM\\SpeciesList\\SpeciesList.csv')
spPyr = unique (spPyr$species)
toEliminate = which(spPyr=='')
spPyr = spPyr [(-1) * toEliminate]

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
mySettings <- occTest::defaultSettings()
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

dir.create("D:/Write/Pep/PyrSDM/OccData")
spPyrDone = list.dirs(path ="D:/Write/Pep/PyrSDM/OccData" ,full.names = F,recursive = F)
spPyrNotDone = setdiff(spPyr,spPyrDone)
cl = makeCluster(5)
registerDoSNOW(cl)
iterations = length(spPyrNotDone)
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
lp = loaded_packages()

out = foreach::foreach(spModel = spPyrNotDone,.options.snow=opts, .errorhandling = 'remove',.packages = lp$package) %dopar% {
  
  outDownload  = try({
    print(which (spPyr == spModel) )
    occRGBIF = rgbif::occ_search(scientificName = spModel,limit = 10000)
    #get the citations for the species
    cit = gbif_citation(occRGBIF)
    #need to go from 2-digit to three digit for occ profiler
    occRGBIF$data$countryCode3c = countrycode::countrycode (occRGBIF$data$countryCode,origin ='iso2c',destination = 'iso3c')
    #write inputs 

    if (user=='serradiaz'){
      jspModel = occProfileR:::.join.spname(spModel)
      pathInputs = paste0('D:/Write/Pep/PyrSDM/OccData/',jspModel)
      }
    if (user=='pserra'){
      pathInputs = paste0('/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/',spModel)
      jspModel = occTest:::.join.spname(spModel)
    }
    dir.create (pathInputs,showWarnings = T,recursive = T)
    save(occRGBIF,file = paste0(pathInputs,'/', jspModel,'.RData'))
    save(cit,file = paste0(pathInputs,'/', jspModel,'_citations.RData'))  }
    ,silent = T)
  
  if (class(outDownload)=="try-error") {return ("Error")} else {return ('Download ended')}
  
}
registerDoSEQ()
stopCluster(cl)
spPyrDone = list.dirs(path ="D:/Write/Pep/PyrSDM/OccData" ,full.names = F,recursive = F)
jspPyr  =  unlist (lapply (spPyr, function (x) occProfileR:::.join.spname(x)))
spPyrNotDone = setdiff(jspPyr,spPyrDone)
idsSpNotDone = which (! jspPyr %in% spPyrDone)
spPyrNotDone = spPyr[idsSpNotDone]





#this species could not be downloaded by rgbif, we use BIEN
spPyrNotDoneRGBIF = spPyrNotDone
lapply (1:length(spPyrNotDoneRGBIF), function (i){
  print (i)
  spModel= spPyrNotDoneRGBIF[i]
  occBIEN = BIEN_occurrence_species(species = spModel,cultivated = T,only.new.world = F,native.status = T,political.boundaries = T,collection.info = T)
  if (nrow (occBIEN)==0) {return (NULL) ; print ('No records found')}
  #need to go from 2-digit to three digit for occ profiler
  occBIEN$country = countrycode::countrycode (occBIEN$country ,origin ='country.name',destination = 'iso3c')
  #write inputs 
  if (user=='serradiaz'){
    jspModel = occProfileR:::.join.spname(spModel)
    pathInputs = paste0('D:/Write/Pep/PyrSDM/OccData/',jspModel)
  }
  if (user=='serradiaz'){
    jspModel = occTest:::.join.spname(spModel)
    pathInputs = paste0('/Users/pserra/RS/cleanOccAnalysis/QAQCinputs/',jspModel)
  }  
  dir.create (pathInputs,showWarnings = F,recursive = T)
  save(occBIEN,file = paste0(pathInputs,'/', jspModel,'_BIEN.RData'))
})

#recheck folder names
pathInputs = "D:/Write/Pep/PyrSDM/OccData"
ldirs = list.dirs(path = pathInputs,full.names = T,recursive = F)
lapply (ldirs, function (x){
  ssplit = strsplit (x, ' ')[[1]]
  if (length(ssplit)) {
    ssplpit2 = paste(ssplit,collapse = '_')
    file.rename(from = x,to = ssplpit2)
  }
})


#check occurrence counts 
ldirs = list.dirs(path = pathInputs,full.names = T,recursive = F)
library (parallelsugar)
NoccSp= mclapply (ldirs,mc.cores= 10,FUN =  function (l){
  print (which (ldirs==l))
  
  allfiles = list.files(path = l,full.names = T)
  dfile = allfiles[1]
  spName = strsplit (basename(dfile),split = '.RData')[[1]]
  occDat = load(dfile)
  occDat = (get(occDat))
  if (class(occDat)=='gbif') {occDat = occDat$data; occinfo='RGBIF'} else {occinfo='BIEN'}
  Nocc=nrow(occDat)
  data.frame(spName, Nocc)
  
})
NoccSp = do.call (rbind, NoccSp)


#how many occurrence per species ?
Nclassification = cut (NoccSp$Nocc,c(10,50,200,500,100001))
table (Nclassification)
spLarge = as.character (NoccSp[which(NoccSp$Nocc>200),'spName'])


#compute QAQC in occProfileR parallel MODE for those species with >200 records
out = lapply (spLarge[2:6], function (spModel){
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
    
    if (nrow (occRGBIF$data) <200) {return(spModel)}

    #LOAD THE DATA 
    pathInputs = "D:/Write/Pep/PyrSDM/OccData"
    allfiles = list.files(path = paste0(pathInputs,'/',spModel) ,full.names = T)
    dfile = allfiles[1]
    spName = strsplit (basename(dfile),split = '.RData')[[1]]
    occDat = load(dfile)
    occDat = (get(occDat))
    if (class(occDat)=='gbif') {occDat = occDat$data; occinfo='RGBIF'} else {occinfo='BIEN'}

    
    #build the Settings for QAQC
    rgbifSettings <- occTest:::defaultSettings()
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
    spProfile_1 = occProfileR:::occurrenceClassify(sp.name= spSci,
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
