### SECONDARY FUNCTIONS TO DETECT OUTLIERS ####
# presPercentile ====
#' Identify percentile of presences
#'
#' @param xy  x y coordinates. Class sf points
#' @param percent numeric. Defaults to 95. Percentage to def object
#' @param unin character vector. 
#' @param unout character vector.
#' @author Josep M Serra Diaz
#' @keywords internal 
#' @description  Divide raster by the sum of all cells.

.presPercentile=function (xy,
                         percent = 95,
                         unin = c("m", "km"),
                         unout = c("ha","km2", "m2")) {

  if (!inherits(xy, "sf") )
    xy <- sf::st_as_sf (xy)
  if (!inherits(xy, "sf") )
    stop("xy should be of class sf")
  pfs <- sf::st_crs(xy)
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

  if (min(table(id)) < 5) stop("At least 5 relocations are required to fit a home range")
  id <- factor(id)
  xy <- as.data.frame(sf::st_coordinates(xy))
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
    xy.t <- st_as_sf(xy.t,coords = c(1,2),crs=4326)
    return(list(xy.t=xy.t,dist.from.centroid=di))
  })
  res
}


# mcpOutliers ====
#' @title Find outlying occurrence data in geographical space from convex hulls
#' @description Computes minimum convex hulls and identifies outlying points
#' @details only works for WGS84 lat long 
#' @param data \emph{data.frame}  with coordinate information
#' @param xcoord \emph{character}. Column with X coordinate information in data
#' @param ycoord \emph{character}. Column with Y coordinate information in data
#' @param percentage \emph{numeric}. Percentage of data used to build the minimum convex polygon
#' @keywords internal
#' @export
#' @examples
#' myPres<-read.csv(system.file('./ext/exampleOccData.csv',package='occTest'))|>
#'   dplyr::select(c('MAPX','MAPY'))
#' mcp_test(data = myPres,xcoord = 'MAPX',ycoord = 'MAPY')
#' @importFrom stats na.omit 
#' @seealso \link[occTest]{geoOutliers}
#' @returns Returns a \emph{logical} vectors detecting outliers
#' @author Coline CF Boonmann <colineboonman@@bio.au.dk>

mcp_test <-  function(data,xcoord,ycoord,percentage=95){
  rownames(data) = seq(1:nrow(data))
  xy = na.omit(data[,c(xcoord,ycoord)])
  xy = xy[!duplicated(xy),]
  if(nrow(xy)>=5){
    
    xy = sf::st_as_sf(xy,coords = c(xcoord,ycoord), crs = sf::st_crs(4326))
    xy$id = "JB_was_here"
    MCP = adehabitatHR::mcp(sf::as_Spatial(xy),percent=percentage)|>
      sf::st_as_sf()
    xy$irow = 1:nrow(xy)
    xy = suppressWarnings(sf::st_intersection(xy,MCP))
    xy = as.data.frame(sf::as_Spatial(xy))
    outliers = c(!(rownames(data) %in% xy$irow))
    return(outliers)
  }else{
    warning("N<5 bservations is not enough to run mcp outlier test.")
    return (rep(NA,nrow (xy)))
    }
}
