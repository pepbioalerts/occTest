#### load some needed packages and functions
library (occTest)
library(raster)
library (sf)

#load data
occ.species.nogeoscrub = system.file('ext/SampleData/Sp3Occurrence_v4.csv'
                                     ,package='occTest')
occ.species.nogeoscrub <-  read.table (file = occ.species.nogeoscrub,header = T,sep = ',',as.is = T)


countries.pol

countries.pol = system.file('ext/CountryLevel',package='occTest')
countries.pol = readOGR(dsn = '/Users/pepserra/Documents/RS/occTest/Rpackage/inst/ext/CountryLevel',layer='Countries_SPDF_MED')


ras.env = system.file('ext/AllEnv.tif',package='occTest')
ras.env = raster::stack (ras.env)

ras.dem = system.file('ext/DEM.tif',package='occTest')
ras.dem = raster (ras.dem)

ras.humaninfluence = system.file('ext/HII.tif',package='occTest')
ras.humaninfluence = raster(ras.humaninfluence)



#prepare projections and plot data
occ.sp <- occ.species.nogeoscrub[,c('MAPX','MAPY')]
occ.sp <- occ.sp [complete.cases(occ.sp),]
coordinates (occ.sp) <-  ~ MAPX + MAPY
occ.sp@proj4string <- sp:::CRS (projection (ras.dem))
plot(ras.dem)
plot(occ.sp, cex=1,add=T)

#inspect data information in occurrence table
head(occ.species.nogeoscrub)
tail(occ.species.nogeoscrub)
###############################################################################
#### PROFILING
###############################################################################

  #### EXAMPLE 1: RUN FUNCTION WITH THE MINIMUM PARAMETERS

  #here we obviate the informationon time, country of record, elevation...just presence and absence
  #note here that we indeed add a native country (Spain) and an invasive country (France)
  #we add all the

occ.profiled.1=occurrenceProfile(output.dir = getwd(),
                                 sp.table = occ.species.nogeoscrub,
                                 sp.name = "my Species",
                                 r.env = ras.env,
                                 r.dem = ras.dem ,
                                 countries.shapefile = NULL,
                                 countryfield.shapefile = 'ISO',
                                 ntv.ctry = c('ESP'),
                                 inv.ctry = c('FRA') ,
                                 ras.hii = ras.humaninfluence,
                                 x.field = 'MAPX',
                                 y.field = 'MAPY',a.field = c("UNCERTAINTY_X_M","UNCERTAINTY_Y_M"),e.field='ELEVATION')


#The result is a list.
str (occ.profiled.1)
#The first element is the full table with all the data
head (occ.profiled.1[[1]])
names (occ.profiled.1[[1]])
#The first element is the short version with only coordinates and quality labels
head (occ.profiled.1[[2]])

#take a look at the diferent quality grades
table(occ.profiled.1$occ_short_profile$quality.grade)
#take a look at the diferent qualifier tags
table(occ.profiled.1$occ_short_profile$qualifiers)
#take a look at the diferent quality lables
table(occ.profiled.1$occ_short_profile$quality.label)

#map the results
full.qaqc <-occ.profiled.1$occ_full_profile
proposed.color.grading <- data.frame (row.names =
                                        c('A','B','C','D','E','F','G','H'),
                                      color.qgrade=c('#4575b4','#74add1',
                                                     '#abd9e9','#e0f3f8',
                                                     '#fee090','#fdae61',
                                                     '#f46d43','#d73027'))
color.to.plot<-as.character(sapply(full.qaqc$quality.grade,
                                   function (x){
                                     proposed.color.grading[x,'color.qgrade'] }))
plot (ras.dem)
points (x =full.qaqc$MAPX,y=full.qaqc$MAPY,col='black',pch=22,bg=color.to.plot)

#distribution of the quality of the plots
library (ggplot2)
dat.to.plot <- full.qaqc
dat.to.plot$grade.color<- color.to.plot
ggplot (aes(x=quality.grade,fill=quality.grade),data = dat.to.plot)+
  geom_bar() +
  ggtitle('Example 1')


#### EXAMPLE 2: RUN FUNCTION ADDING THE INFORMATION IN THE DATABASE ON ELEVATION, TIME, AND LOCALITY

occ.profiled.2=occurrenceProfile(output.dir =
                                   "/Users/pepserra/Documents/RS/occTest_sampleData/",
                                 sp.table = occ.species.nogeoscrub,
                                 sp.name = "My species",
                                 r.env = ras.env,
                                 r.dem = ras.dem ,
                                 countries.shapefile = countries.pol,
                                 countryfield.shapefile = 'ISO',
                                 ntv.ctry = c('ESP'),
                                 inv.ctry = c('FRA') ,
                                 ras.hii = ras.humaninfluence,
                                 x.field = 'MAPX' ,
                                 y.field = 'MAPY',
                                 #added info
                                 t.field = 'DATE',
                                 e.field = 'ELEVATION',
                                 l.field = 'LOCALITYNAME')




#we know have added one more quality grade

# (because now we can track botanical gardens (e.g. new class E together with missing environment))
table(occ.profiled.1$occ_short_profile$quality.grade)
table(occ.profiled.2$occ_short_profile$quality.grade)

#we know also have more quality tags added
# (we have the tag T, meaning we have a timestamp in the record, and the tag e, meaining that the elevation recorded and the elevation extracted is less than 100m (default))
table(occ.profiled.1$occ_short_profile$qualifiers)
table(occ.profiled.2$occ_short_profile$qualifiers)

#take a look at the diferent quality lables
table(occ.profiled.2$occ_short_profile$quality.label)

#map the results
  #full.qaqc.2 <-occ.profiled.2$occ_full_profile
  # color.to.plot.2 <- as.character  (sapply (full.qaqc.2$quality.grade, function (x){proposed.color.grading[x,'color.qgrade'] }))
  # plot (ras.dem)
  # points (x =full.qaqc.2$MAPX ,y=full.qaqc.2$MAPY ,col='black',pch=22,cex=0.1,bg=color.to.plot,text(full.qaqc.2$MAPY~full.qaqc.2$MAPX, labels =full.qaqc.2$quality.grade,pos=2))

#distribution of the quality of the plots
library (ggplot2)
dat.to.plot <- full.qaqc
dat.to.plot$grade.color<- color.to.plot
ggplot (aes(x=quality.grade,fill=quality.grade),data = dat.to.plot)+geom_bar() +ggtitle('Example 2')


#### EXAMPLE 3: RUN FUNCTION ADDING THE INFORMATION OF COUNTRY RECORD

#adding the information on country record collected may dramatically change the results
#the grading system will hugely downgrade any recored with a wrong or ABSENT country of reference name.
#we recommend only to enable this feature if all the countries of the records are provided, as any NA will be considred as a 'very low quality' record

occ.profiled.3    <- occurrenceProfile  (output.dir ="/Users/pepserra/Documents/RS/occTest_sampleData/",
                                         sp.table = occ.species.nogeoscrub,
                                         sp.name = "My species",
                                         r.env = ras.env,
                                         r.dem = ras.dem ,
                                         countries.shapefile = countries.pol,
                                         countryfield.shapefile = 'ISO',
                                         ntv.ctry = c('ESP'),
                                         inv.ctry = c('FRA') ,
                                         ras.hii = ras.humaninfluence,
                                         x.field = 'MAPX' ,
                                         y.field = 'MAPY',
                                         t.field = 'DATE',e.field = 'ELEVATION',l.field = 'LOCALITYNAME',
                                         #we now added info on the country record
                                         c.field = "COUNTRYRECORD"
)

#check: only three classes now, almost all corresponding to category F
full.qaqc <-occ.profiled.3$occ_full_profile
ggplot (aes(x=quality.grade,fill=quality.grade),data = full.qaqc)+geom_bar() +ggtitle('Example 2')


#### EXAMPLE 4: RUN FUNCTION without information on native range
#setting to NULL the ntv.ctry will make the grading system shift, as now everything is considered native range
occ.profiled.4    <- occurrenceProfile   (output.dir = "/Users/pepserra/Documents/RS/occTest_sampleData/",
                                         sp.table = occ.species.nogeoscrub,
                                         sp.name = "My species",
                                         r.env = ras.env,
                                         r.dem = ras.dem ,
                                         countries.shapefile = countries.pol,
                                         countryfield.shapefile = 'ISO',
                                         #changed info on ntv.ctry
                                         ntv.ctry = NULL ,
                                         inv.ctry = c('FRA') ,
                                         ras.hii = ras.humaninfluence,
                                         x.field = 'MAPX' ,
                                         y.field = 'MAPY',
                                         t.field = 'DATE',e.field = 'ELEVATION',l.field = 'LOCALITYNAME'

)
full.qaqc.4 <-occ.profiled.4$occ_full_profile

#check: only three classes now, almost all corresponding to category F
#note that the point from Andorra has incrased
#also note, however, that this changes the way the range is computed, as well as the potential outliers.
#for instance, the inclusion of the point of Andorra changed also the categorization of geographical outliers, because we include this point now in the alphahull calculation
#therefore, accepting all countries as native countries can have impacts on the grading system.

#check
ex2 <- full.qaqc.2[,c('ID_GRAFIC','quality.label')]
names (ex2) <- c('ID_GRAFIC','Q.2')
ex4 <- full.qaqc.4[,c('ID_GRAFIC','quality.label')]
names (ex4) <- c('ID_GRAFIC','Q.4')
plyr::join (ex2,ex4)

#### EXAMPLE 5: RUN FUNCTION without information invasive range

#setting to NULL THE ntv.ctry will make the grading system to classify potential ranges
#it can thus affect outlier detection, as those records are not part of the outlier detection.
occ.profiled.5    <- occurrenceProfile (output.dir = "/Users/pepserra/Documents/RS/occTest_sampleData/",
                                         sp.table = occ.species.nogeoscrub,
                                         sp.name = "My species",
                                         r.env = ras.env,
                                         r.dem = ras.dem ,
                                         countries.shapefile = countries.pol,
                                         countryfield.shapefile = 'ISO',
                                         #changed info on ntv.ctry
                                         ntv.ctry = c('ESP') ,
                                         inv.ctry = NULL ,
                                         ras.hii = ras.humaninfluence,
                                         x.field = 'MAPX' ,
                                         y.field = 'MAPY',
                                         t.field = 'DATE',e.field = 'ELEVATION',l.field = 'LOCALITYNAME',excludeUnknownRanges = F,excludeNotmatchCountry = F

)

#check: only three classes now, almost all corresponding to category F
#note that the point from Andorra has incrased
#also note, however, that this changes the way the range is computed, as well as the potential outliers.
#for instance, the inclusion of the point of Andorra changed also the categorization of geographical outliers, because we include this point now in the alphahull calculation
#therefore, accepting all countries as native countries can have impacts on the grading system.

#check, a lot of F appeard (unknown range)
full.qaqc.5 <-  occ.profiled.5[[1]]
ex5 <- full.qaqc.5[,c('ID_GRAFIC','quality.label')]
names (ex5) <- c('ID_GRAFIC','Q.5')
#analyze de difference between the two
plyr::join (ex4,ex5)




