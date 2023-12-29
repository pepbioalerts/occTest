# landSeaFilter ====
#' @title Filter occurrence by habitat (terrestrial/non-Terrestrial)
#' @description Filter the occurrence records according to whether they should be in land masses or not
#' @param df \emph{data.frame} of species occurrences
#' @param xf \emph{character}. The field in df containing the x coordinates
#' @param yf \emph{character}. The field in df containing the y coordinates
#' @param habType \emph{character}. Define the species habitat. Only "terrestrial" and "sea" implemented. 
#' @param habPol \emph{polygon sf}. Multipolygon object indicating land masses
#' @param verbose \emph{logical} Print messages? Default TRUE
#' @returns \emph{list} of filtered and non-filtered occurrrence recprds by land use type
#' @details The user can define the habitat type and provide a an input landmass polygon
#' @keywords filter
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples \dontrun{
#' xyDF <- data.frame (x=c(0,42),y=c(0,1),Reason=NA)
#' landSeaFilter(df=xyDF,xf='x',yf='y')
#' }
#' @export

landSeaFilter             =function(df,
                                    xf,
                                    yf,
                                    habType=NULL,verbose=TRUE,habPol=NULL) {

  
  #load high-res land masses
  if (is.null(habPol)){
    dest_url <- 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/allLand10.rds'
    outFile = paste0(tempdir(),'/allLand10.rds')
    if (!file.exists(outFile)) utils::download.file(url=dest_url,destfile = outFile)
    land = readRDS (outFile)
  } else { land = habPol }
  #select coordinates and load sf points
  xydat <- df[,c(xf,yf)]
  pts = sf::st_as_sf(xydat,coords=c(1,2),crs=terra::crs(land))
  #intersect
  intersectMatrix <- sf::st_intersects(x = pts, y = land,sparse = TRUE)
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
