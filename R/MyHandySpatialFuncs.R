### Spatial functions


#FUNCTION TO GET THE ADEQUATE PROJECTION AND THEN CREATE A BUFFER
#An AEQ projection centred on each point will project equal distances in all directions, such as a buffered circle. 
#(A Mercator projection will have some distortions in north and eastern directions, since it wraps around the side of a cylinder.)
#from : https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world

aeqd.buffer <- function(p, r)
{
  stopifnot(length(p) == 1)
  aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  p@coords[[2]], p@coords[[1]])
  projected <- spTransform(p, CRS(aeqd))
  buffered <- gBuffer(projected, width=r, byid=TRUE)
  spTransform(buffered, p@proj4string)
}



bbox2extent <- function (obj){
  stopifnot(class(obj) == 'matrix')
  
  extent.object <- extent(obj['x','min'],
                          obj['x','max'],
                          obj['y','min'],
                          obj['y','max'])
  
  return (extent.object)
  
  
}





#FUNCTION TO GET FROM ALHPA HULLS TO SHAPEFILES

ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
  require(alphahull)
  require(maptools)
  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  # Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf,xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0){
    # Convert each arc to a line segment
    linesj <- list()
    prevx<-NULL
    prevy<-NULL
    j<-1
    for(i in 1:nrow(xdf)){
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)){
        prevx<-x
        prevy<-y
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
        if (i == nrow(xdf)){
          #We have got to the end of the dataset
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
          prevx[length(prevx)]<-prevx[1]
          prevy[length(prevy)]<-prevy[1]
          coordsj<-cbind(prevx,prevy)
          colnames(coordsj)<-NULL
          # Build as Line and then Lines class
          linej <- Line(coordsj)
          linesj[[j]] <- Lines(linej, ID = as.character(j))
        } else {
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
        }
      } else {
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)]<-prevx[1]
        prevy[length(prevy)]<-prevy[1]
        coordsj<-cbind(prevx,prevy)
        colnames(coordsj)<-NULL
        # Build as Line and then Lines class
        linej <- Line(coordsj)
        linesj[[j]] <- Lines(linej, ID = as.character(j))
        j<-j+1
        prevx<-NULL
        prevy<-NULL
      }
    }
    # Promote to SpatialLines
    lspl <- SpatialLines(linesj)
    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}
