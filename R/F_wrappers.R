#### FUNCTIONAL WRAPPERS
###### ========================================

#'Runs tests and validates data
#'
#'
#' @param spOcc data.frame. Object with the coordinate data.
#' @param speciesName character. Name of the species.
#' @param env raster or rasterStack. Environmental data (e.g. typically climatic).
#' @param x name of the field with the coordinate x. Default 'x'
#' @param y name of the field with the coordinate y. Default 'y'
#' @param data name of the field with the values of the reported timestamp of the record. Default NULL
#' @param isoCountry name of the field with the values of the reported country of the record. Default NULL
#' @param classification character. Indicates the thresholds phylosofphy applied to classify errors in occurrence data. Possible values 'strict','relaxed','custom'
#' @return a list of two. First element is a dataframe with profiled occurrence records with their associated profiled labels. Second element is a dataframe with all outputs of the analysis implemented.
#' @note #There are several parameters in the function. The majority of them can be adjusted, but we also provide default values. We recommend those default values if the user is to use the geospatial data included in the package.
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export

#wallace simple wrapping funciton
occSimpleClassification = function (spOcc,env,speciesName='My species',x='x',y='y',date=NULL,isoCountry=NULL,classification=NULL ){
  
  #set up params
  mySettings = occProfileR::defaultSettings()
  mySettings$tableSettings$x.field <- x
  mySettings$tableSettings$y.field <- y
  mySettings$tableSettings$t.field <- date
  mySettings$tableSettings$c.field <- isoCountry
  
  
  #run test functions 
  output =  occurrenceTests   (sp.name = speciesName,
                               sp.table = spOcc, 
                               r.env = env,
                               tableSettings =mySettings$tableSettings)
  
  #classify
  if(classification == 'custom') print ('not implemented yet')
  
  
  output 
  
  
}
