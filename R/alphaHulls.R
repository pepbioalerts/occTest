### Alpha Hull Analysis related scripts

# getPointsOutAlphaHull ====
#' @title Generate Spatial object from Alpha Hull object
#' @description  Generate polygon based on alpha hulls at a given alpha parmater 
#' @details 
#' Based on rangeBuilder::getDynamicRange without cropping to sea and without increments
#' \code{alpha = initialAlpha}, and will then increase \code{alpha} by
#' \code{alphaIncrement} until both the \code{fraction} and \code{partCount}
#' conditions are met.
#' 
#' If the conditions cannot be satisfied, then a minimum convex hull is
#' returned.
#' 
#' If \code{clipToCoast} is set to "terrestrial" or "aquatic", the resulting
#' polygon is clipped to the coastline, using the dataset
#' provided with this package.
#' 
#' @param x dataframe of coordinates in decimal degrees, with a minimum of 3
#' rows.
#' @param alpha the starting value for alpha.
#' @param coordHeaders the column names for the longitude and latitude
#' columns, respectively.  If x has two columns, these are assumed to be
#' longitude and latitude, and \code{coordHeaders} is ignored.
#' @param proj the projection information for x. The default is currently the
#' only supported option.
#' @param alphaCap Max alpha value before function aborts and returns a
#' minimum convex hull.
#' @param verbose logical. print messages? Default to FALSE
#' @return a list with 2 elements: \item{hull}{ a SpatialPolygons object }
#' \item{alpha}{ the alpha value that was found to satisfy the criteria.  If a
#' convex hull was returned, this will list MCH.  }
#' @author Pascal Title (original version), Josep M Serra-Diaz (modifs)
#' @seealso Alpha hulls are created with \code{\link{ahull}}.
#' @examples
#' @importFrom methods slot<-
#' @keywords internal

getPointsOutAlphaHull <- function(x,  alpha = 2, coordHeaders = c('Longitude', 'Latitude'), 
                                  #buff = 1000, parameter not implemented
                                  proj = "+proj=longlat +datum=WGS84",  
                                  verbose = FALSE, alphaCap = 20) {
  
  if (proj != "+proj=longlat +datum=WGS84") {
    stop("Currently, proj can only be '+proj=longlat +datum=WGS84'.")
  }
  
  if (ncol(x) == 2) {
    coordHeaders <- c(1,2)
  }
  
  #reduce to unique coordinates [this is probably done before, but it is a safety measure]
  x <- x[!duplicated(x[,coordHeaders]), coordHeaders]
  x <- x[stats::complete.cases(x),]
  
  if (nrow(x) < 3) {
    stop('This function requires a minimum of 3 unique coordinates (after removal of duplicates).')
  }
  
  #Alpha hulls cannot be generated if first 3 points are linear. 
  while ((x[1,coordHeaders[1]] == x[2,coordHeaders[1]] & x[2, coordHeaders[1]] == x[3, coordHeaders[1]]) | (x[1,2] == x[2, coordHeaders[2]] & x[2, coordHeaders[2]] == x[3, coordHeaders[2]])) {
    x <- x[sample(1:nrow(x),size = nrow(x)),]
  }
  
  #create spatialpoints
  x <- sp::SpatialPoints(x, proj4string = sp::CRS(proj))
  x <- sp::remove.duplicates(x)
  if (length(x) < 3) {
    stop('This function requires a minimum of 3 unique coordinates.')
  }
  
  #create alpha hull and calculate fraction of occurrences that fall within
  #continue until fraction is reached
  problem <- FALSE
  if (verbose) {cat('\talpha:', alpha, '\n')}
  
  #check alpha requested makes sense
  if (alpha > alphaCap) {
    stop (paste0('Alpha parmeter',alpha,'is bigger than Alpha Cap',alphaCap))
  }
  
  #perfomr hull
  hull <- try(alphahull::ahull(data.frame(x),alpha = alpha), silent = TRUE)
  if (inherits(hull, 'try-error')) {
    stop('Alpha hull not built')
  }
  
  hull <- try ( .ah2sp(hull, proj4string = sp::CRS('+proj=longlat +datum=WGS84')), silent=TRUE)
  if (!is.null(hull)) {
      slot(hull, "polygons") <- lapply(slot(hull, "polygons"),  .checkPolygonsGEOS2)
  }
  
  if (is.null(hull) | inherits(hull, 'try-error') | !cleangeo::clgeo_IsValid(hull)) {
    stop ('Not valid GEOS for alphaHull built ')
    }
  
  #how many points are within hull?
  slot(hull, "polygons") <- lapply(slot(hull, "polygons"),  .checkPolygonsGEOS2)
  if (!cleangeo::clgeo_IsValid(hull)) stop ('Invalid GEOS for alphahull built')
  
  pointWithin <- rgeos::gIntersects(x, hull, byid = TRUE)
  pointWithin <- pointWithin[1,]
  pointsOut <- !pointWithin * 1
  return (pointsOut)
  
}



### .ah2sp  ======
#' @title Convert Alpha Hull object into a shapefile
#' @details Function written by Andrew Bevan, found on R-sig-Geo, and modified by Pascal Title
#' @param x an alpha hull object
#' @param increment numeric. Increments
#' @param rnd numeric. Decimal rounding
#' @param proj4string crs object with the spatial projection.
#' @param tol numeric. tolerance
#' @return a sp polygon object
#' @author Pascal Title (original version), Josep M Serra-Diaz (modifications)
#' @seealso Alpha hulls are created with \code{\link{ahull}}.
#' @importFrom methods slot<-
#' @importFrom sp Lines Polygons

.ah2sp <- function(x, increment=360, rnd=10, proj4string=sp::CRS(as.character(NA)),tol=1e-4) {
  if (!inherits(x, "ahull")) {
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  
  #correct for possible arc order strangeness (Pascal Title addition 29 Nov 2013)
  k <- 1
  xdf <- cbind(xdf, flip = rep(FALSE, nrow(xdf)))
  repeat{
    if (is.na(xdf[k+1, 'end1'])) {
      break
    }
    #cat(k, '\n')
    if (xdf[k,'end2'] == xdf[k+1,'end1']) {
      #cat('case 1\n')
      k <- k + 1
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & !xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & !xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 2\n')
      k <- k + 1
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & !xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 3\n')
      m <- which(xdf$end1[k+1:nrow(xdf)] == xdf[k,'end2']) + k
      xdf <- rbind(xdf[1:k,],xdf[m,],xdf[setdiff((k+1):nrow(xdf),m),])
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & !xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 4\n')
      m <- which(xdf$end2[k+1:nrow(xdf)] == xdf[k,'end2']) + k
      tmp1 <- xdf[m,'end1']
      tmp2 <- xdf[m,'end2']
      xdf[m,'end1'] <- tmp2
      xdf[m,'end2'] <- tmp1
      xdf[m,'flip'] <- TRUE
      xdf <- rbind(xdf[1:k,], xdf[m,], xdf[setdiff((k+1):nrow(xdf), m),])
    } else {
      k <- k + 1
    }
  }	
  
  
  # Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf, xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0) {
    # Convert each arc to a line segment
    linesj <- list()
    prevx <- NULL
    prevy <- NULL
    j <- 1
    for(i in 1:nrow(xdf)) {
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2), 0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- alphahull::anglesArc(v, theta)
      if (rowi['flip'] == TRUE){ angles <- rev(angles) }
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)) {
        prevx <- x
        prevy <- y
        # added numerical precision fix (Pascal Title Dec 9 2013)
      } else if ((x[1] == round(prevx[length(prevx)],rnd) | abs(x[1] - prevx[length(prevx)]) < tol) && (y[1] == round(prevy[length(prevy)],rnd) | abs(y[1] - prevy[length(prevy)]) < tol)) {
        if (i == nrow(xdf)){
          #We have got to the end of the dataset
          prevx <- append(prevx ,x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
          prevx[length(prevx)] <- prevx[1]
          prevy[length(prevy)] <- prevy[1]
          coordsj <- cbind(prevx,prevy)
          colnames(coordsj) <- NULL
          # Build as Line and then Lines class
          linej <- sp::Line(coordsj)
          linesj[[j]] <- sp::Lines(linej, ID = as.character(j))
        } else {
          prevx <- append(prevx, x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
        }
      } else {
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)] <- prevx[1]
        prevy[length(prevy)] <- prevy[1]
        coordsj <- cbind(prevx,prevy)
        colnames(coordsj)<-NULL
        # Build as Line and then Lines class
        linej <- sp::Line(coordsj)
        linesj[[j]] <- Lines(linej, ID = as.character(j))
        j <- j + 1
        prevx <- NULL
        prevy <- NULL
      }
    }
    
    #Drop lines that will not produce adequate polygons (Pascal Title addition 9 Dec 2013)
    badLines <- vector()
    for (i in 1:length(linesj)){
      if (nrow(linesj[[i]]@Lines[[1]]@coords) < 4){
        badLines <- c(badLines,i)
      }
    }
    if (length(badLines) > 0){linesj <- linesj[-badLines]}
    
    # Promote to SpatialLines
    lspl <- sp::SpatialLines(linesj)
    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start and end points that are the same
    lns <- methods::slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- methods::slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- sp::SpatialPolygons(list(sp::Polygons(lapply(list_of_Lines, function(x) { sp::Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- sp::SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}

# .checkPolygonsGEOS2 ====
#' @title  Check polygon geometry
#' @description
#' @details inspired provided by  maptools package and from P Title in rangeBuilder
#' @param obj an alpha hull object
#' @param properly logic.
#' @param force logic.
#' @param useSTRtree logic.
#' @return a sp polygon object
#' @author Pascal Title (original version), Josep M Serra-Diaz (modifications)
#' @seealso Alpha hulls are created with \code{\link{ahull}}. \cr
#' see maptools and RangeBuilder package

.checkPolygonsGEOS2 <- function(obj, properly = TRUE, force = TRUE, useSTRtree = FALSE) {
  if (!is(obj, "Polygons")) 
    stop("not a Polygons object")
  comm <- try(rgeos::createPolygonsComment(obj), silent = TRUE)
  if (!inherits(comm, "try-error") && !force) {
    comment(obj) <- comm
    return(obj)
  }
  pls <- slot(obj, "Polygons")
  IDs <- slot(obj, "ID")
  n <- length(pls)
  if (n < 1) 
    stop("Polygon list of zero length")
  uniqs <- rep(TRUE, n)
  if (n > 1) {
    if (useSTRtree) 
      tree1 <- rgeos::gUnarySTRtreeQuery(obj)
    SP <- sp::SpatialPolygons(lapply(1:n, function(i) Polygons(list(pls[[i]]), ID = i)))
    for (i in 1:(n - 1)) {
      if (useSTRtree) {
        if (!is.null(tree1[[i]])) {
          res <- try(rgeos::gEquals(SP[i, ], SP[tree1[[i]],], byid = TRUE), silent = TRUE)
          if (inherits(res, "try-error")) {
            warning("Polygons object ", IDs, ", Polygon ", i, ": ", res)
            next
          }
          if (any(res)) {
            uniqs[as.integer(rownames(res)[res])] <- FALSE
          }
        }
      }
      else {
        res <- try(rgeos::gEquals(SP[i, ], SP[uniqs,], byid = TRUE), silent = TRUE)
        if (inherits(res, "try-error")) {
          warning("Polygons object ", IDs, ", Polygon ", i, ": ", res)
          next
        }
        res[i] <- FALSE
        if (any(res)) {
          wres <- which(res)
          uniqs[wres[wres > i]] <- FALSE
        }
      }
    }
  }
  if (any(!uniqs)) 
    warning(paste("Duplicate Polygon objects dropped:", paste(wres, collapse = " ")))
  pls <- pls[uniqs]
  n <- length(pls)
  if (n < 1) 
    stop("Polygon list of zero length")
  if (n == 1) {
    oobj <- Polygons(pls, ID = IDs)
    comment(oobj) <- rgeos::createPolygonsComment(oobj)
    return(oobj)
  }
  areas <- sapply(pls, slot, "area")
  pls <- pls[order(areas, decreasing = TRUE)]
  oholes <- sapply(pls, function(x) slot(x, "hole"))
  holes <- rep(FALSE, n)
  SP <- sp::SpatialPolygons(lapply(1:n, function(i) Polygons(list(pls[[i]]), ID = i)))
  if (useSTRtree) 
    tree2 <- rgeos::gUnarySTRtreeQuery(SP)
  for (i in 1:(n - 1)) {
    if (useSTRtree) {
      if (!is.null(tree2[[i]])) {
        if (properly) 
          res <- rgeos::gContainsProperly(SP[i, ], SP[tree2[[i]], ], byid = TRUE)
        else res <- rgeos::gContains(SP[i, ], SP[tree2[[i]], ], byid = TRUE)
      }
      else {
        res <- FALSE
      }
    }
    else {
      if (properly) 
        res <- rgeos::gContainsProperly(SP[i, ], SP[-(1:i), ], byid = TRUE)
      else res <- rgeos::gContains(SP[i, ], SP[-(1:i), ], byid = TRUE)
    }
    wres <- which(res)
    if (length(wres) > 0L) {
      nres <- as.integer(rownames(res))
      holes[nres[wres]] <- !holes[nres[wres]]
    }
  }
  for (i in 1:n) {
    if (oholes[i] != holes[i]) 
      pls[[i]] <- sp::Polygon(slot(pls[[i]], "coords"), hole = holes[i])
  }
  oobj <- Polygons(pls, ID = IDs)
  comment(oobj) <- rgeos::createPolygonsComment(oobj)
  oobj
}


