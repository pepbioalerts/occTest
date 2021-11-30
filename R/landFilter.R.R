#' @title Filter occurrence by habitat (terrestrial/non-Terrestrial)
#'
#' @description Filter the occurrence recoreds  according to whether they should be in land masses or not
#' @details
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @return list
#' @keywords filter
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export

landSeaFilter             =function(df,
                                    xf,
                                    yf,
                                    geom= NULL,
                                    habType=NULL,verbose=T) {

  
  # df=dat; xf=x.field; yf=y.field; .ntv.ctry=ntv.ctry;
  #select coordinates and load sf points
  xydat <- df[,c(xf,yf)]
  pts = sf::st_as_sf(xydat,coords=c(1,2),crs=crs(land))
  #load high-res land masses
  land = readRDS (system.file('ext/land/allLand10.rds',package='occTest'))
  #intersect
  intersectMatrix <- sf::st_intersects(x = pts, y = land,sparse = T)
  intersectMatrix <- Matrix::as.matrix(intersectMatrix)
  ptsInLand <- Matrix::rowSums(intersectMatrix)
  #points in out of the land/sea
  if (is.null(habType)) {
    habType = 'terrestrial'
    if (verbose) message ('Assuming a terrestrial species. Modify parameter Habitat otherwise')}
  if (habType=='terrestrial') idRows <- !as.logical(ptsInLand)
  if (habType=='sea')         idRows <- as.logical(ptsInLand)
  #output
  df$Exclude <- idRows * 1
  df$Reason [idRows] <- paste('not',habType,'habitat',collapse = ' ')
  df$coordIssues_wrongHabitat_value = idRows * 1
  df$coordIssues_wrongHabitat_test = idRows
  out <- list (stay = df[which(df$Exclude==1),], continue = df[which(df$Exclude!=1),])
  return (out)
  
}