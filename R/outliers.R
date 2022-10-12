
# presPercentile ====
#' Identify percentile of presences
#'
#' @param xy  x y coordinates. Class SpatialPoints
#' @param percent numeric. Defaults to 95. Percentage to def object
#' @param unin character vector. 
#' @param unout character vector.
#' @author Josep M Serra Diaz
#' @keywords internal 
#' @description  Divide raster by the sum of all cells.

presPercentile=function (xy,
                         percent = 95,
                         unin = c("m", "km"),
                         unout = c("ha","km2", "m2")) {

  if (!inherits(xy, "SpatialPoints"))
    stop("xy should be of class SpatialPoints")
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
      warning("xy should contain only one column (the id of the animals), id ignored")
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    } else {
      id <- xy[[1]]
    }
  } else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }

  if (min(table(id)) < 5) stop("At least 5 relocations are required to fit an home range")
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
      acons <- key[di <= quantile(di, percent/100)]
    } else { acons=key }
    xy.t <- df.t[acons, ]
    sp::coordinates(xy.t)=c(1,2)
    return(list(xy.t=xy.t,dist.from.centroid=di))
  })
  res
}

# findSpatialOutliers =====
#' Find outlying occurrence data in geographical space
#'
#' @param myPres  raster* object
#' @param pvalSet numeric. p value set to identify outliers
#' @param checkPairs logical.  
#' @param verbose logical. print messages
#' @author Cory Merow
#' @keywords internal
#' @description  Divide raster by the sum of all cells.
#' @return a numeric vector indicating which rows are spatial outliers.
#' @examples  
#' k <- data.frame (x=c(runif (n = 100),1000),y=c(runif (n = 100),1000))
#' k <- sp::SpatialPoints(k)
#' occTest::findSpatialOutliers(k)
#' @export
findSpatialOutliers=function(myPres,
                             pvalSet=1e-5,
                             checkPairs=TRUE,verbose=TRUE){

  pres.inliers=myPres
  sp.toss.coord=NULL
  pval=0
  #toss singles
  while(pval<pvalSet){
    dists=presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    gt=outliers::grubbs.test(dists)
    #gt=dixon.test(dists)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    (pval=gt$p.value)
    # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
    if(gt$p.value<pvalSet){
      toss=which.max(dists)
      # IDs in the original data frame
      sp.toss.coord=rbind(sp.toss.coord,sp::coordinates(pres.inliers)[toss,])
      pres.inliers=pres.inliers[-toss,]
      dists=dists[-toss]
    }
  }
  # toss pairs
  if(checkPairs){
    if(length(pres.inliers)<31){
      pval=0
      # By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
      #while(pval<pvalSet){
      dists=presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
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

  if(!is.null(sp.toss.coord)){
    coor=sp::coordinates(myPres)
    sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
  } else {sp.toss.id=NULL}
  if (verbose)
    print(paste0(length(sp.toss.id),' geographic outliers found'))
  sp.toss.id
}

# findEnvOutliers =====
#' @title Find outlying occurrence data in environmental space
#' @description Environmental outliers
#' @details
#' @param myPres a `SpatialPointsDataFrame`
#' @param myEnv a `RasterStack` of env layers. If NULL, it is assumed that `myPres` is a data.frame of the environmental values (columns) at presence locations (rows)
#' @param checkPairs logical. Default to FALSE (TRUE not implemented).
#' @param pvalSet numeric; p-value used in Grubb's test for outlier (see package `outliers`)
#' @param verbose logic. Should messages be printed out?
#' @keywords internal
#' @export
#' @return Returns a list of SpatialPointsDataFrames with (1) good presence points (2) spatial outliers and (3) environmental outliers.
#' @author Cory Merow <cory.merow@@gmail.com>

findEnvOutliers=function(myPres,
                         myEnv=NULL,
                         pvalSet=1e-5,
                         checkPairs=FALSE,verbose=TRUE){
  #  for testing
  #  myPres=presDF; pvalSet=1e-5; checkPairs=FALSE; myEnv=NULL
  #  myEnv=env
  #  myPres=myPres; env=myEnv; pvalSet=1e-5
  if(!is.null(myEnv)){ p.env=raster::extract(myEnv,myPres)} else {p.env=myPres}
  p.env=scale(p.env)
  # remove variables that are the same for all observations
  f=which(apply(p.env,2,function(x) !all(is.nan(x))))
  p.env=p.env[,f]
  pres.inliers=p.env
  if (any (class (p.env) == 'numeric')) row.id = as.character(p.env)
  if (any (class (p.env) %in% c('matrix','tibble','data.frame')))    row.id=apply( p.env, 1 , paste , collapse = "-" )
  env.toss.id=NULL
  pval=0
  while(pval<pvalSet){
    if (any (class (p.env) == 'numeric')) dists = p.env
    if (any (class (p.env) %in% c('matrix','tibble','data.frame'))) dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    gt=outliers::grubbs.test(dists)
    pval=gt$p.value
    # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
    if(gt$p.value<pvalSet){
      toss=which.max(dists)
      
    
      # IDs in the original data frame
      if (!any (class (p.env) == 'numeric')) {
        thisID=paste(p.env[toss,],collapse='-')
        env.toss.id=c(env.toss.id,which(row.id == thisID))
        p.env=p.env[-toss,]  
      }
      if (any (class (p.env) == 'numeric')) {
        thisID=as.character(p.env[toss])
        env.toss.id=c(env.toss.id,which(row.id == thisID))
        p.env=p.env[-toss]  
      }
    
    }
  }
  if(checkPairs) print('checkPairs not yet implemented')
  env.toss.id=unique(env.toss.id)
  if (verbose)
    print(paste0(length(env.toss.id),' environmental outliers found'))
  env.toss.id
}
