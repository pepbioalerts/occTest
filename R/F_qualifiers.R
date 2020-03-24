#### FUNCTION QUALIFIERS (FOR TAGGING)


### qualifiers of range type
.qualifier.invasive.range <- function  (df,
                                       tag='i',
                                       qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){
  #first scope
  if (sum (df$quality.grade %in% qualifier.lab.scp)==0) {out <- rep (NA,nrow(df)); return (out)}
  
  #determining
  range.info <- df$countryStatusRangeAnalysis_wrongINV_test
  range.info [is.na(range.info)] = 1
  range.info = ifelse (range.info == 0,tag,NA)
  
  #scoping
  range.info[! df$quality.grade %in% qualifier.lab.scp]=NA
  
  return (range.info)
}

.qualifier.native.range <- function  (df,
                                     tag='n',
                                     qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){
  #first scope
  if (sum (df$quality.grade %in% qualifier.lab.scp)==0) {out <- rep (NA,nrow(df)); return (out)}
  
  #determining
  range.info <- df$countryStatusRangeAnalysis_wrongNTV_test
  range.info [is.na(range.info)] = 1
  range.info = ifelse (range.info == 0,tag,NA)
  
  #scoping
  range.info[! df$quality.grade %in% qualifier.lab.scp]=NA
  
  return (range.info)
}

#qualifiers of timestamp
.qualifier.timestamp <- function  (df,
                                  tf=t.field ,
                                  tag='t',
                                  qualifier.lab.scp = c('A','B','C','D','E'), verbose=F){

  #first scope
  if (sum (df$quality.grade %in% qualifier.lab.scp)==0) {out <- rep (NA,nrow(df)); return (out)}
  
  # #determining
  timestamp <- ifelse (is.na(df[,tf]) | (as.character (df[,tf]) == '') , NA, tag)
  
  #scoping
  timestamp[! df$quality.grade %in% qualifier.lab.scp]=NA
  
  return (timestamp)

}

#qualifiers of precision
.qualifier.precision <- function (df,
                                 tag='p',
                                 xf=x.field,
                                 yf=y.field,
                                 .s=res.in.minutes,
                                 .e=res.in.minutes,
                                 qualifier.lab.scp = c('A','B','C','D','E'), verbose=F ) {

  #first scope
  if (sum (df$quality.grade %in% qualifier.lab.scp)==0) {out <- rep (NA,nrow(df)); return (out)}
  
  #determining
  datpc <- biogeo::precisioncheck (dat= df, x = xf,y = yf,s=.s, e=.e)
  preci <- ifelse (datpc$preci == 1,NA,tag)
  
  #scoping
  preci[! df$quality.grade %in% qualifier.lab.scp]=NA
  
  return (preci)



}

#qualifier of elevation threshold
.qualifier.elevation <- function (df=dat,
                                 qualifier.lab.scp = c('A','B','C','D','E'), verbose=F,
                                 tag='e'
                                 ) {

  #first scope
  if (sum (df$quality.grade %in% qualifier.lab.scp)==0) {out <- rep (NA,nrow(df)); return (out)}
  
  all.NA.test = length(is.na(df$elevDiff_test)) == nrow (df) 
  name.not.exist = ! ('elevDiff_test' %in% names(df))
  if (any(all.NA.test,name.not.exist)) {out = rep (NA,length.out=nrow(df)) ;  return(out)}
  
  #remerging
  out <- ifelse (join.df$df$elevDiff_test == 1, tag, NA)
  #scoping
  out[! df$quality.grade %in% qualifier.lab.scp]=NA
  
  
  return (out)


}


