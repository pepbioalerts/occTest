#example with tables from species examples from Wallace

library(spocc)
df <- occ(query = 'Martes martes')
occ.data <-occ2df(df)

#we change the names here but we could just change it in the tableSettings
names (occ.data)[2] <- 'x'
names (occ.data)[3] <- 'y'


#load needed data
library (raster)
rrr = raster::getData(name='worldclim',var='bio', res=10,path = "D:/Tmp")

#classify occurrences
library(occProfileR)
out = occurrenceClassify(sp.name='Martes martes',sp.table = occ.data,r.env = rrr,interactiveMode = T)
out = occurrenceClassify(sp.table = occ.data,r.env = rrr,interactiveMode = F)

occurrenceClassify()

head(out$occ_short_profile)
table(out$occ_short_profile$quality.grade)

