#### occOutliers env.r =====

#' @title Find outlying occurrence data in environmental space
#' @description Environmental outliers
#' @details
#' Perform outlier detaction in the environmental space for interquartile range (boxplot method), 
#' Dixon (see \link[outliers]{dixon.test}), Rosner test (see Dixon (see \link[EnvStats]{rosnerTest})),
#' and Grubbs (see Dixon (see \link[outliers]{grubbs.test})) tests. 
#' @param pres an \emph{sf} points object describing the locations of species records with their associated environmenal data. 
#' @param method \emph{character}. Methods to identify outlying points options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet \emph{numeric}. User-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs \emph{logical}. check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner \emph{numeric}. Number between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
#' @returns a vector with record indices considered outliers. 
#' @export
#' @examples
#' myPres=read.csv(system.file('ext/exampleOccData.csv',
#'                             package='occTest'))
#' #select coordinates only
#' myPres <- myPres[,2:3]
#' myPres <- myPres[complete.cases(myPres),]
#' #extract values
#' myEnv <- terra::rast(system.file('ext/AllEnv.tif',package='occTest'))
#' myPres <- sf::st_as_sf(myPres,coords = c(1,2))
#' myPresDF <- cbind (myPres,data.frame(terra::extract(myEnv,myPres)))
#' #find outliers
#' presOut=findEnvOutliers(pres=myPresDF,pvalSet=1e-5)
#' @return Returns the indices of environmental outliers
#' @author Cory Merow <cory.merow@@gmail.com>
#' @importFrom EnvStats rosnerTest
#' @seealso \link[occTest]{envOutliers} \link[occTest]{findOutlyingPoints} \link[outliers]{dixon.test} \link[outliers]{grubbs.test} \link[EnvStats]{rosnerTest}
#' @references 
#' Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.
#' Dixon, W.J. (1951). Ratios involving extreme values. Ann. Math. Stat. 22, 1, 68-78.
#' Rorabacher, D.B. (1991). Statistical Treatment for Rejection of Deviant Values: Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95 percent Confidence Level. Anal. Chem. 83, 2, 139-146.
#' Grubbs, F.E. (1950). Sample Criteria for testing outlying observations. Ann. Math. Stat. 21, 1, 27-58.
#' Rosner, B. (1975). On the Detection of Many Outliers. Technometrics 17, 221–227.
#' Rosner, B. (1983). Percentage Points for a Generalized ESD Many-Outlier Procedure. Technometrics 25, 165–172.
findEnvOutliers=function(pres,
                         #myEnv=NULL,
                         pvalSet=1e-5,
                         method='grubbs',
                         checkPairs=TRUE,
                         kRosner=NULL){
  p.env=sf::st_drop_geometry(pres)
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
      if(nrow(p.env)<3) {break} 
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
        pval=gt$p.value
        # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
        if(gt$p.value<pvalSet){
          toss=utils::tail(order(dists),2)
          env.toss.id = toss
          # IDs in the original data frame
          # sp.toss.coord=rbind(sp.toss.coord, sp::coordinates(pres.inliers)[toss,])
          # pres.inliers=pres.inliers[-toss,]
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
    if(kRosner >= length(dists)) {
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
#' @description Spatial outliers
#' @details
#' Perform outlier detaction in the environmental space for interquartile range (boxplot method), 
#' Dixon (see \link[outliers]{dixon.test}), Rosner test (see Dixon (see \link[EnvStats]{rosnerTest})),
#' and Grubbs (see Dixon (see \link[outliers]{grubbs.test})) tests.
#' @param pres an \emph{sf} points object describing the locations of species records with their associated environmenal data. 
#' @param method \emph{character}. Methods to identify outlying points options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet \emph{numeric}. User-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs \emph{logical}. check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner \emph{numeric}. Number between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
#' @returns a vector with record indices considered outliers. 
#' @export
#' @examples
#' myPres=read.csv(system.file('./ext/exampleOccData.csv',
#'                             package='occTest')) |> dplyr::select(c('MAPX','MAPY'))
#' myPres = sf::st_as_sf(myPres[complete.cases(myPres),],
#'                       coords =c('MAPX','MAPY'))
#' findSpatialOutliers(pres=myPres, pvalSet=1e-5)
#' 
#' @return Returns the indices of environmental outliers
#' @author Cory Merow <cory.merow@@gmail.com>
#' @importFrom EnvStats rosnerTest
#' @seealso \link[occTest]{geoOutliers} \link[occTest]{findOutlyingPoints} \link[outliers]{dixon.test} \link[outliers]{grubbs.test} \link[EnvStats]{rosnerTest}
#' @references 
#' Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.
#' Dixon, W.J. (1951). Ratios involving extreme values. Ann. Math. Stat. 22, 1, 68-78.
#' Rorabacher, D.B. (1991). Statistical Treatment for Rejection of Deviant Values: Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95 percent Confidence Level. Anal. Chem. 83, 2, 139-146.
#' Grubbs, F.E. (1950). Sample Criteria for testing outlying observations. Ann. Math. Stat. 21, 1, 27-58.
#' Rosner, B. (1975). On the Detection of Many Outliers. Technometrics 17, 221–227.
#' Rosner, B. (1983). Percentage Points for a Generalized ESD Many-Outlier Procedure. Technometrics 25, 165–172.

findSpatialOutliers=function(pres,
                             pvalSet=1e-5,
                             method='grubbs',
                             checkPairs=FALSE,
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
    while(pval<pvalSet & nrow(pres.inliers)>3){
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
        sp.toss.coord=rbind(sp.toss.coord,sf::st_coordinates(pres.inliers)[toss,])
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
          #stop()
        }
        if(gt$p.value<pvalSet){
          toss=utils::tail(order(dists),2)
          # IDs in the original data frame
          sp.toss.coord=rbind(sp.toss.coord,sf::st_coordinates(pres.inliers)[toss,])
          pres.inliers=pres.inliers[-toss,]
        }
        #}
      }	
    }
    if(!is.null(sp.toss.coord)){
      coor=sf::st_coordinates(pres)
      sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
    } 
  } # end grubbs
  
  if(any(method=='iqr')) {
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    sp.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
    if(nrow(pres.inliers)<3 | nrow(pres.inliers)>30 ){
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
    if(kRosner >= length(dists)) {
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
#' @description This function is a wrapper to identify spatial or environmental outliers using different methods
#' @details
#' Perform outlier detaction in the environmental space for interquartile range (boxplot method), 
#' Dixon (see \link[outliers]{dixon.test}), Rosner test (see Dixon (see \link[EnvStats]{rosnerTest})),
#' and Grubbs (see Dixon (see \link[outliers]{grubbs.test})) tests.
#' @param pres an \emph{sf} points object describing the locations of species records with their associated environmenal data. 
#' @param spOutliers logical; perform spatial outlier analysis
#' @param envOutliers logical; perform environmental outlier analysis
#' @param method \emph{character}. Methods to identify outlying points options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pval \emph{numeric}. User-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs \emph{logical}. check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner \emph{numeric}. Number between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
#' @param verbose \emph{logical}; print messages?
#' @export
#' @examples
#' #select coordinates only
#' myPres=read.csv(system.file('ext/exampleOccData.csv',
#'                             package='occTest'))
#' #select coordinates only
#' myPres <- myPres[,2:3]
#' myPres <- myPres[complete.cases(myPres),]
#' #extract values
#' myEnv <- terra::rast(system.file('ext/AllEnv.tif',package='occTest'))
#' myPres <- sf::st_as_sf(myPres,coords = c(1,2))
#' myPresDF <- cbind (myPres,data.frame(terra::extract(myEnv,myPres)))
#' #find outliers
#' findOutlyingPoints(pres=myPresDF,pval=1e-5)
#' @returns Returns an \emph{sf} polygon object  with additional logical columns indicating spatial outliers (`spOutliers`) and environmental outliers ('`envOutliers`).
#' @author Cory Merow <cory.merow@@gmail.com>
#' @seealso \link[occTest]{geoOutliers} \link[occTest]{findEnvOutliers} \link[occTest]{findSpatialOutliers} \link[outliers]{dixon.test} \link[outliers]{grubbs.test} \link[EnvStats]{rosnerTest}
#' @references 
#' Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.
#' Dixon, W.J. (1951). Ratios involving extreme values. Ann. Math. Stat. 22, 1, 68-78.
#' Rorabacher, D.B. (1991). Statistical Treatment for Rejection of Deviant Values: Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95 percent Confidence Level. Anal. Chem. 83, 2, 139-146.
#' Grubbs, F.E. (1950). Sample Criteria for testing outlying observations. Ann. Math. Stat. 21, 1, 27-58.
#' Rosner, B. (1975). On the Detection of Many Outliers. Technometrics 17, 221–227.
#' Rosner, B. (1983). Percentage Points for a Generalized ESD Many-Outlier Procedure. Technometrics 25, 165–172.
findOutlyingPoints=function(pres,
                            spOutliers=TRUE,
                            envOutliers=TRUE,
                            method='grubbs',
                            pval=1e-5,
                            checkPairs=FALSE,
                            kRosner=5,
                            verbose=TRUE){

  if(spOutliers) { 
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
.vgrep=function(pattern,x){mapply(function(y){grep(y,pattern)},x)}
