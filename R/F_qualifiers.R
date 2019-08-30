#### FUNCTION QUALIFIERS (FOR TAGGING)


### qualifiers of range type
.qualifier.invasive.range <- function  (df,
                                       tag='i',
                                       qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){
  #scoping
  df2 <- subset (df,subset= df$quality.grade %in% qualifier.lab.scp)
  if (nrow (df2) == 0) {out <- rep (NA,nrow(df)); return (out)}
  #determining
  range.info <- df2$countryStatusRangeAnalysis_wrongINV_test
  range.info [is.na(range.info)] <- 1
  range.info <- ifelse (range.info == 0,tag,NA)
  df2$range.info <- range.info
  #remerging
  join.df <-  suppressMessages(plyr:::join(df,df2))
  out <- join.df$range.info
  return (out)
}

.qualifier.native.range <- function  (df,
                                     tag='n',
                                     qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){
  #scoping
  df2 <- subset (df,subset= df$quality.grade %in% qualifier.lab.scp)
  if (nrow (df2) == 0) {out <- rep (NA,nrow(df)); return (out)}
  #determining
  range.info <- df2$countryStatusRangeAnalysis_wrongNTV_test
  range.info [is.na(range.info)] <- 1
  range.info <- ifelse (range.info == 0,tag,NA)
  df2$range.info <- range.info
  #remerging
  join.df <- suppressMessages( plyr:::join(df,df2))
  out <- join.df$range.info
  return (out)

}

#qualifiers of timestamp
.qualifier.timestamp <- function  (df,
                                  tf=t.field ,
                                  tag='t',
                                  qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){

  #scoping
  df2 <- subset (df,subset= df$quality.grade %in% qualifier.lab.scp)
  if (nrow (df2) == 0) {out <- rep (NA,nrow(df)); return (out)}
  #determining
  df2$timestamp <- ifelse (is.na(df2[,tf]) | (df2[,tf] == '') , NA, tag)
  #remerging
  join.df <-  suppressMessages(plyr:::join(df,df2))
  out <- join.df$timestamp
  return (out)

}

#qualifiers of precision
.qualifier.precision <- function (df,
                                 tag='p',
                                 xf=x.field,
                                 yf=y.field,
                                 .s=res.in.minutes,
                                 .e=res.in.minutes,
                                 qualifier.lab.scp = c('A','B','C','D','E'), verbose=F ) {

  #scoping
  df2 <- subset (df,subset= df$quality.grade %in% qualifier.lab.scp)
  if (nrow (df2) == 0) {out <- rep (NA,nrow(df)); return (out)}
  #determining
  datpc <- biogeo::precisioncheck (dat= df2, x = xf,y = yf,s=.s, e=.e)
  df2$preci <- ifelse (datpc$preci == 1,NA,tag)
  #remerging
  join.df <-  suppressMessages( plyr:::join(df,df2))
  out <- join.df$preci
  return (out)



}

#qualifier of elevation threshold
.qualifier.elevation <- function (df=dat,
                                 qualifier.lab.scp = c('A','B','C','D','E'), verbose=F,
                                 tag='e'
                                 ) {

  
  #determining
  df2 <- subset (df,subset= df$quality.grade %in% qualifier.lab.scp)
  
  all.NA.test = length(is.na(df2$elevDiff_test)) == nrow (df2) 
  name.not.exist = ! ('elevDiff_test' %in% names(df2))
  if (any(all.NA.test,name.not.exist)) {out = rep (NA,length.out=nrow(df)) ;  return(out)}
  
  #remerging
  df2$elevdiff2 <- df2$elevDiff_test
  join.df <-  suppressMessages(plyr:::join(df,df2))
  out <- ifelse (join.df$elevdiff2 == 1, tag, NA)
  return (out)


}


