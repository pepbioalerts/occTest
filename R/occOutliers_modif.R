#### occOutliers env.r =====

#' @title Find outlying occurrence data in environmental space
#'
#' @description Environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' myEnv=raster::stack(system.file('extdata/AllEnv.tif',package='occOutliers'))
#' names(myEnv)=read.table(system.file('extdata/layerNames.csv',package='occOutliers'))[,1]
#' myPresDF=sp::SpatialPointsDataFrame(myPres,data.frame(raster::extract(myEnv,myPres)))
#' presOut=findEnvOutliers(pres=myPresDF,pvalSet=1e-5)
#' @return Returns the indices of environmental outliers
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findEnvOutliers=function(pres,
                         #myEnv=NULL,
                         pvalSet=1e-5,
                         method='grubbs',
                         checkPairs=TRUE,
                         kRosner=NULL){
  #  for testing
  #  pres=presDF; pvalSet=1e-5; checkPairs=F; myEnv=NULL
  #  myEnv=env
  #  env=myEnv; pvalSet=1e-5
  #if(!is.null(myEnv)){ p.env=raster::extract(myEnv,pres)
  #} else {p.env=pres}
  p.env=pres@data
  p.env=base::scale(p.env)
  sp.toss.coord=NULL
  sp.toss.id=NULL
  # remove variables that are the same for all observations
  f=which(apply(p.env,2,function(x) !all(is.nan(x))))
  p.env=p.env[,f,drop=F]
  pres.inliers=p.env
  row.id=apply( p.env, 1 , paste , collapse = "-" )
  env.toss.id=NULL
  #dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
  if(any(method=='grubbs')){
    pval=0
    #tmp.dists=dists
    while(pval<pvalSet & length(pres.inliers)>3){
      # do this within the loop to recompute the centroid ofter each outlier is removed
      dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
      # this is close enough to fixed, in the outer function it'll report that everything was tossed
      if(nrow(p.env)<3) break 
      gt=outliers::grubbs.test(dists)
      pval=gt$p.value
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        thisID=paste(p.env[toss,],collapse='-')
        env.toss.id=c(env.toss.id,which(row.id == thisID))
        #tmp.dists=tmp.dists[-toss]
        p.env=p.env[-toss,,drop=F]
      }  
    }
    
    #if(checkPairs) print('checkPairs not yet implemented')
    # toss pairs
    if(checkPairs){ # taken from spatial - must test
      if(length(pres.inliers)<31 & length(pres.inliers)>3){
        pval=0
        dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
        # By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
        #while(pval<pvalSet){
        # duplicated
        #dists=.presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
        gt=outliers::grubbs.test(dists,type=20)
        
        #gt=dixon.test(dists)
        #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
        (pval=gt$p.value)
        # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
        if(gt$p.value<pvalSet){
          toss=utils::tail(order(dists),2)
          # IDs in the original data frame
          sp.toss.coord=rbind(sp.toss.coord, sp::coordinates(pres.inliers)[toss,])
          pres.inliers=pres.inliers[-toss,]
        }
        #}
      }	
    }
  }# end grubbs
  
  if(any(method=='iqr')){
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    env.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
    if(length(pres.inliers)<3 | length(pres.inliers)>30 ){
      warning('dixon test only applies to sample sizes between [3,30]. skipping this taxon')
      return(env.toss.id)
    }
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    if(length(unique(dists))==1) {
      warning('all records are the same distance from the environmental centroid so outliers cannot be detected. maybe your records come from gridded data. skipping this taxon')
      return(env.toss.id)
    }
    dt=outliers::dixon.test(dists,type=0,two.sided = FALSE)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    if(dt$p.value<pvalSet) env.toss.id=which.max(dists)
  }
  
  if(any(method=='rosner')){
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    if(kRosner <= length(dists)) {
      warning('kRosner must be an integer less than the number of presence records, skipping this taxon')
      return(env.toss.id)
    }
    rt=EnvStats::rosnerTest(dists,kRosner,alpha=pvalSet)
    if(any(rt$all.stats$Outlier)) env.toss.id=utils::tail(order(dists),kRosner)
  }
  
  # if(!is.null(sp.toss.coord)){ is this needed from spatial?
  #   coor=coordinates(pres)
  #   sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
  # } else {sp.toss.id=NULL}
  
  unique(env.toss.id)
}



#### occOutliers spatial.r ====
#' @title Find outlying occurrence data in geographic space
#'
#' @description Spatial outliers
#'
#' @details
#' See Examples.
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' presOut=findSpatialOutliers(pres=myPres, pvalSet=1e-5)
#' 
#' @return Returns the indices of spatial outliers
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findSpatialOutliers=function(pres,
                             pvalSet=1e-5,
                             method='grubbs',
                             checkPairs=TRUE,
                             kRosner=NULL){
  #  for testing
  #   pvalSet=1e-5; checkPairs=T
  pres.inliers=pres
  sp.toss.coord=NULL
  sp.toss.id=NULL
  if(any(method=='grubbs')){
    pval=0
    #tmp.dists=dists
    #toss singles
    while(pval<pvalSet & length(pres.inliers)>3){
      # i want to recompute the distances once an outlier is removed because the outlier biased the centroid of the group, which influences distances
      dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
      gt=outliers::grubbs.test(dists)
      #gt=dixon.test(dists)
      #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
      (pval=gt$p.value)
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      if(is.na(pval)) {
        warning('p-value for grubbs test was NA. sample size is too small to implement the test')
        break
      }
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        sp.toss.coord=rbind(sp.toss.coord,sp::coordinates(pres.inliers)[toss,])
        pres.inliers=pres.inliers[-toss,]
        #tmp.dists=tmp.dists[-toss]
      }
    }
    if(length(dists)<3) warning('All but two records were deemed outliers. The Grubbs test may not be appropriate for these data.')
    # toss pairs
    if(checkPairs){
      if(length(pres.inliers)<31 & length(pres.inliers)>3){
        pval=0
        #tmp.dists=dists
        dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
        
        # By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
        #while(pval<pvalSet){
        gt=outliers::grubbs.test(dists,type=20)
        
        #gt=dixon.test(dists)
        #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
        (pval=gt$p.value)
        # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
        if(is.na(pval)) {
          warning('p-value for grubbs test checking for pairs of outliers was NA. sample size is too small to implement the test')
          break
        }
        if(gt$p.value<pvalSet){
          toss=utils::tail(order(dists),2)
          # IDs in the original data frame
          sp.toss.coord=rbind(sp.toss.coord, sp::coordinates(pres.inliers)[toss,])
          pres.inliers=pres.inliers[-toss,]
        }
        #}
      }	
    }
    if(!is.null(sp.toss.coord)){
      coor=sp::coordinates(pres)
      sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
    } 
  } # end grubbs
  
  if(any(method=='iqr')) {
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    sp.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
    if(length(pres.inliers)<3 | length(pres.inliers)>30 ){
      warning('dixon test only applies to sample sizes between [3,30]. skipping this taxon')
      return(sp.toss.id)
    }
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(length(unique(dists))==1) {
      warning('all records are the same distance from the spatial centroid so outliers cannot be detected. maybe your records come from gridded data. skipping this taxon')
      return(sp.toss.id)
    }
    dt=outliers::dixon.test(dists,type=0,two.sided = FALSE)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    if(dt$p.value<pvalSet) sp.toss.id=which.max(dists)
  }
  
  if(any(method=='rosner')){
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(kRosner <= length(dists)) {
      warning('kRosner must be an integer less than the number of presence records, skipping this taxon')
      return(sp.toss.id)
    }
    rt=EnvStats::rosnerTest(dists,kRosner,alpha=pvalSet)
    if(any(rt$all.stats$Outlier)) sp.toss.id=utils::tail(order(dists),kRosner)
  }
  
  sp.toss.id
}
#### occOutliers findOutliers.r =====
#' @title Find outlying occurrence data
#'
#' @description Spatial or environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param spOutliers logical; perform spatial outlier analysis
#' @param envOutliers logical; perform environmental outlier analysis
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'. Only a single value can be used; if mutliple tests are desired, call `findOuliers` multiple times specifying different methods.
#' @param pval user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
#' @param verbose logical; print messages
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' myEnv=raster::stack(system.file('extdata/AllEnv.tif',package='occOutliers'))
#' names(myEnv)=read.table(system.file('extdata/layerNames.csv',package='occOutliers'))[,1]
#' myPresDF=sp::SpatialPointsDataFrame(myPres,data.frame(raster::extract(myEnv,myPres)))
#' presOut=findOutlyingPoints(pres=myPresDF,
#'                            spOutliers=TRUE,
#'                            envOutliers=TRUE,
#'                            pval=1e-5)
#' 
#' @return Returns the SpatialPointsDataFrame provided with additional logical columns indicating spatial outliers (`spOutliers`) and environmental outliers ('`envOutliers`).
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


findOutlyingPoints=function(pres,
                            spOutliers=TRUE,
                            envOutliers=TRUE,
                            method='grubbs',
                            pval=1e-5,
                            checkPairs=TRUE,
                            kRosner=3,
                            verbose=TRUE){
  
  #  for testing
  #  pres=myPresDF;  verbose=T
  #  envOutliers=T; spOutliers=T;  pval=1e-5; 
  if(!any(class(pres)==c('SpatialPoints','SpatialPointsDataFrame'))) stop('Please make your presence data a SpatialPoints or SpatialPointsDataFrame object and try again')
  if(spOutliers) { 
    
      if (class(pres)==c('SpatialPoints')) {
        #get them to become spdfs
        pres  = sp::SpatialPointsDataFrame(pres,
                                           data.frame(spOutlier= rep(NA,length(pres)))
        )
        
      }
    sp.toss.id=findSpatialOutliers(pres=pres,pvalSet=pval,checkPairs=checkPairs,
                                   method=method,kRosner=kRosner)
    
    if(is.null(sp.toss.id)) pres$spOutlier=NA
    if(!is.null(sp.toss.id)) {
      pres$spOutlier=FALSE
      pres$spOutlier[sp.toss.id]=TRUE
    }
    
    if(verbose) print(paste0(length(sp.toss.id),' geographic outlier(s) found with method ',method))
  } else {sp.toss.id=NULL}
  
  if(envOutliers) {
    if (class(pres)==c('SpatialPoints')) {
      #get them to become spdfs
      pres  = sp::SpatialPointsDataFrame(pres,
                                         data.frame(envOutlier= rep(NA,length(pres)))
      )
      
    }
    env.toss.id=findEnvOutliers(pres=pres,pvalSet=pval,checkPairs=checkPairs,
                                method=method,kRosner=kRosner)
    if(is.null(env.toss.id)) pres$envOutlier=NA
    if(!is.null(env.toss.id)) {
      pres$envOutlier=FALSE
      pres$envOutlier[env.toss.id]=TRUE
    }
    if(verbose) print(paste0(length(env.toss.id),' environmental outlier(s) found with method ',method))
    if(nrow(pres)-length(env.toss.id) < 2){
      warning(paste0('pretty much all the presences came up as environmental outliers. The only case Ive seen this in was when there were two obvious outliers and all the other presences had the same exact environment. In any case, I cant find any outliers.'))
    }
  }
  
  return(pres)
}



#### occOutliers







#### occOutliers utils.r =====
#' omit outlying pres
#' @param xy data.frame with 2 columns
#' @param percent numeric on [0,100]
#' @param unin character; input units. option are "m" and "km"
#' @param unout character; input units. option are "ha","km2", "m2"
# @param ... additional functions to be passed to \code{\link[raster]{writeRaster}}
# @description  Divide raster by the sum of all cells.
# @export
# unin = c("m");  unout ='m2'
.presPercentile=function (xy, percent = 95, unin = c("m", "km"), unout = c("ha","km2", "m2")) {
  # for testing
  # xy=pres; percent = NULL; unin = c("m", "km"); unout = c("ha","km2", "m2")
  #if (!inherits(xy, "SpatialPoints")) # not needed; i already check it
  #stop("xy should be of class SpatialPoints")
  if (ncol(sp::coordinates(xy)) > 2)
    stop("xy should be defined in two dimensions")
  pfs <- sp::proj4string(xy)
  if(!is.null(percent)){
    if (length(percent) > 1)
      stop("only one value is required for percent")
    if (percent > 100) {
      warning("Using all relocations (percent>100)")
      percent <- 100
    }
  }
  unin <- match.arg(unin)
  unout <- match.arg(unout)
  if (inherits(xy, "SpatialPointsDataFrame")) {
    if (ncol(xy) != 1) {
      #warning("xy should contain only one column (the id of the animals), id ignored")
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    } else {
      id <- xy[[1]]
      if (all(is.na(id))) id <- factor(rep("a", length(id)))
    }
  } else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }
  
  if (min(table(id)) < 4) stop("must have 4 records to proceed")
  id <- factor(id)
  xy <- as.data.frame(sp::coordinates(xy))
  r <- split(xy, id)
  est.cdg <- function(xy) apply(xy, 2, mean)
  cdg <- lapply(r, est.cdg)
  levid <- levels(id)
  res =lapply(1:length(r), function(i) {
    k <- levid[i]
    df.t <- r[[levid[i]]]
    cdg.t <- cdg[[levid[i]]]
    dist.cdg <- function(xyt) {
      d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
      return(d)
    }
    di <- apply(df.t, 1, dist.cdg)
    key <- c(1:length(di))
    if(!is.null(percent)){
      acons <- key[di <= stats::quantile(di, percent/100)]
    } else { acons=key }
    xy.t <- df.t[acons, ]
    sp::coordinates(xy.t)=c(1,2)
    return(list(xy.t=xy.t,dist.from.centroid=di))
  })
  res
}

#' find records mor than 1.5 times the interquartile range beyond the upper quartile 
#' @param dists a numeric vector
.iqrOutlier=function(dists){
  
  #q1=stats::quantile(dists, .25)
  q3=stats::quantile(dists, .25)
  iqr=stats::IQR(dists)
  which(dists > (q3 + 1.5*iqr))
  #subset(df, df$A> (Q1 - 1.5*IQR) & df$A< (Q3 + 1.5*IQR))
}

#' Vectorized version of grep
#' @description vectorized version of grep
#' @param pattern character string containing a regular expression (or character string for ‘fixed = TRUE’) to be matched in the given character vector.  Coerced by ‘as.character’ to a character string if possible.  If a character vector of length 2 or more is supplied, the first element is used with  a warning.
#' @param x a character vector where matches are sought, or an object which can be coerced by ‘as.character’ to a character vector.
# @export
.vgrep=function(pattern,x){mapply(function(y){grep(y,pattern)},x)}
