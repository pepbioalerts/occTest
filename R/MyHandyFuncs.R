#my handy functions

#####
############################################################
join.spname <- function (x){
  a <- strsplit(x,split=' ')
  a <- unlist (a)
  a <- paste (a,collapse = '_')
  return(a)
}

collapse.name <- function (x){
  a <- strsplit(x,split=' ')
  a <- unlist (a)
  a <- paste (a,collapse = '_')
  return(a)
}


#####
############################################################
multiple.strsplit <- function (x,multiple.splits) {
  
  stopifnot( length(multiple.splits)>=1 )
  
  a <- unlist(strsplit (x, split = multiple.splits[1]))
  
  if (length(multiple.splits>1)){
    
    for (i in 2:length(multiple.splits)){
      if (length(a)==1 & a=='') {next()}
      if (is.null(a)) {a <-''; next()}
      a <- unlist( strsplit (a, split = multiple.splits[i] ) )
    }
    
  }
  
  return(a)
}

#####
############################################################

rm.all <- function () {rm (list =  setdiff(ls(),'rm.all') ) ; gc(verbose = F) }


#####
############################################################

lazylogic <- function (e){
  
  o <- eval(parse(text=e))
  ifelse(is.na(o),F,o)
  
}

#####
############################################################

subsetlist.nonNULL <- function (x) {Filter(Negate(is.null), x)}
subsetlist.isNA <- function (x) {Filter(Negate(is.na), x)}


#####
############################################################
unlist.lapply <- function (x,fun) {unlist (lapply (X = x,FUN = fun))}

####
#####################################################################
table.to.dataframe <- function (x){
  n <- names (x)
  nums <- as.numeric (x)
  l <- list()
  for (i in 1:length (n)){
    l[[n[i]]] <- nums[i]
  }
  df <- as.data.frame (l)
  names (df) <- n
  return (df)
}

table.df <- function (x){
  initial.table <- table (x)
  out.df <- as.data.frame (t (as.matrix (initial.table)))
  return (out.df)
}


#####
############################################################
gimme.genus.species <- function (x) {
  a <- strsplit (x, split = ' ')[[1]]
  data.frame (genus=a[1], sp.epith=a[2])
}


### gimme a new version of a file
#####################################

file.new.version <- function (ini.filename, extension, newversion){
  
  if (!is.character(ini.filename) | !is.character(extension) | !is.character(newversion)) {stop('All inputs must be character type')}
  
  new.ext <- paste0('.',extension)
  a <- (strsplit (ini.filename,new.ext)[[1]])
  b <- paste0(a,newversion,new.ext)
  return(b)
}
