#example with tables from species examples from Wallace

library(spocc)
df <- occ(query = 'Martes martes')
occ.data <-occ2df(df)

#we change the names here but we could just change it in the tableSettings
names (occ.data)[2] <- 'x'
names (occ.data)[3] <- 'y'
#we add some countries
occ.data$ctry = coords2country(occ.data[,2:3])
#we make something not work
occ.data$ctry[c(5,55)]<- 'ESP'

#load meta
#load needed data
library (raster)
rrr = raster::getData(name='worldclim',var='bio', res=10,path = "~/RS/tmpData/")


#modify the settings tot take the country into account
mySettings <- defaultSettings()
mySettings$tableSettings$c.field <- 'ctry'
mySettings$tableSettings$t.field <- 'date'


#classify occurrences
library(occTest)

out = occTest(sp.name='Martes martes',
                         sp.table = occ.data,
                         r.env = rrr,
                         interactiveMode = F,
                         tableSettings = mySettings$tableSettings)



nrow (occ.data)
nrow (out$occ_short_profile)
occ.data[1:3,1:3]
View (out$occ_short_profile)


out = occurrenceClassify(sp.table = occ.data,r.env = rrr,interactiveMode = F)


head(out$occ_short_profile)
head(out$occ_full_profile)

table(out$occ_short_profile$quality.grade)

