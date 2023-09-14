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
  
  #create spatial points object
  x <- sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 4326)
  if (nrow(x) < 3) {
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
  
  #perform hull
  hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
  while (inherits(hull, 'try-error') & any(grepl('duplicate points', hull))) {
    ptDist <- sf::st_distance(x, x)
    diag(ptDist) <- NA
    units(ptDist) <- NULL
    closest <- which(ptDist == min(ptDist, na.rm = TRUE), arr.ind = TRUE)
    hull <- try(alphahull::ahull(sf::st_coordinates(x)[- closest[1,1]], alpha = alpha), silent = TRUE)
    }
  if (inherits(hull, 'try-error')) {
    stop('Alpha hull not built')
  }
  hull <- try ( .ah2sf(hull), silent=TRUE)
  validityCheck <- function(hull) {
    if (!is.null(hull) & !inherits(hull, 'try-error')) {
      if (all(sf::st_is_valid(hull))) {
        TRUE
      } else {
        FALSE
      }
    } else {
      FALSE
    }
  }
  if (!all(sf::st_is_valid(hull))) {hull <- sf::st_make_valid(hull)}
  if (! validityCheck(hull)){
    stop('Alpha hull not built')
  }
  
  #how many points are within hull?
  pointWithin <- sf::st_intersects(x, hull,sparse=F)
  if (dim(pointWithin)[2]>1) 
    pointWithin <- as.logical (rowSums (pointWithin))
  
  pointsOut <- as.vector (!pointWithin)*1 
  return (pointsOut)
  
}



### .ah2sf  ======
# Function written by Andrew Bevan, found on R-sig-Geo, and modified by Pascal Title
# modified to support sf objects 17 Nov 2022

.ah2sf <- function(x, increment = 360, rnd = 10, crs = 4326, tol = 1e-4) {
  if (!inherits(x, "ahull")) {
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  
  #correct for possible arc order strangeness (Pascal Title addition 29 Nov 2013)
  k <- 1
  xdf <- cbind(xdf, flip = rep(FALSE, nrow(xdf)))
  repeat {
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
    # linesj <- list()
    linesj <- sf::st_sf(id = 1:nrow(xdf), geometry = sf::st_sfc(lapply(1:nrow(xdf), function(x) sf::st_multilinestring())), crs = 4326)
    prevx <- NULL
    prevy <- NULL
    j <- 1
    for (i in 1:nrow(xdf)) {
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2), 0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- alphahull::anglesArc(v, theta)
      if (rowi['flip'] == TRUE) angles <- rev(angles)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang), rnd)
      y <- round(cc[2] + r * sin(seqang), rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)) {
        prevx <- x
        prevy <- y
        # added numerical precision fix (Pascal Title Dec 9 2013)
      } else if ((x[1] == round(prevx[length(prevx)],rnd) | abs(x[1] - prevx[length(prevx)]) < tol) && (y[1] == round(prevy[length(prevy)],rnd) | abs(y[1] - prevy[length(prevy)]) < tol)) {
        if (i == nrow(xdf)) {
          #We have got to the end of the dataset
          prevx <- append(prevx ,x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
          prevx[length(prevx)] <- prevx[1]
          prevy[length(prevy)] <- prevy[1]
          coordsj <- cbind(prevx,prevy)
          colnames(coordsj) <- NULL
          # Build as Line and then Lines class
          # linej <- Line(coordsj)
          # linesj[[j]] <- Lines(linej, ID = as.character(j))
          linej <- sf::st_linestring(coordsj)
          linesj$geometry[j] <- linej
          
        } else {
          prevx <- append(prevx, x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
        }
      } else {
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)] <- prevx[1]
        prevy[length(prevy)] <- prevy[1]
        coordsj <- cbind(prevx,prevy)
        colnames(coordsj) <- NULL
        # Build as Line and then Lines class
        # linej <- Line(coordsj)
        # linesj[[j]] <- Lines(linej, ID = as.character(j))
        linej <- sf::st_linestring(coordsj)
        linesj$geometry[j] <- linej
        j <- j + 1
        prevx <- NULL
        prevy <- NULL
      }
    }
    
    # drop empty line geometries
    linesj <- linesj[which(sf::st_is_empty(linesj) == FALSE),]
    
    # Drop lines that will not produce adequate polygons (Pascal Title addition 9 Dec 2013, updated 24 May 2023 for sf)
    badLines <- integer()
    for (i in 1:nrow(linesj)) {
      tmp <- sf::st_cast(linesj[i,], 'POINT', warn = FALSE)
      if (nrow(tmp) < 4) {
        badLines <- c(badLines, i)
      }
    }
    
    if (length(badLines) > 0) {
      linesj <- linesj[-badLines, ] 
    }
    
    res <- sf::st_geometry(sf::st_cast(linesj, 'POLYGON'))
    
    
    # # Promote to SpatialLines
    # lspl <- SpatialLines(linesj)
    # # Convert lines to polygons
    # # Pull out Lines slot and check which lines have start and end points that are the same
    # lns <- slot(lspl, "lines")
    # polys <- sapply(lns, function(x) { 
    # crds <- slot(slot(x, "Lines")[[1]], "coords")
    # identical(crds[1, ], crds[nrow(crds), ])
    # }) 
    # # Select those that do and convert to SpatialPolygons
    # polyssl <- lspl[polys]
    # list_of_Lines <- slot(polyssl, "lines")
    # sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), crs=crs)
    # # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    # hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    # areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    # df <- data.frame(hid,areas)
    # names(df) <- c("HID","Area")
    # rownames(df) <- df$HID
    # res <- SpatialPolygonsDataFrame(sppolys, data=df)
    # res <- res[which(res@data$Area > 0),]
  }	
  return(res)
}
