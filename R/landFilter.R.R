#' @title Filter occurrence by habitat (terrestrial/non-Terrestrial)
#'
#' @description Filter the occurrence recoreds  according to whether they should be in land masses or not
#' @details
#' @param df Data.frame of species occurrences
#' @param xf the field in the dataframe containing the x cordinates
#' @param yf the field in the dataframe containing the y cordinates
#' @param habType character. Define the species habitat. Only "terrestrial" and "sea" implented.
#' @param verbose logical. Print messages? Default T
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
                                    habType=NULL,verbose=T) {

  
  #load high-res land masses
  dest_url = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/allLand10.rds'
  outFile = paste0(tempdir(),'/allLand10.rds')
  if (!file.exists(outFile)) utils::download.file(url=dest_url_hii,destfile = outFile)
  land = readRDS (outFile)
  #select coordinates and load sf points
  xydat <- df[,c(xf,yf)]
  pts = sf::st_as_sf(xydat,coords=c(1,2),crs=crs(land))
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
