#### FUNCTIONAL WRAPPERS
###### ========================================
#'Runs tests and validates data
#' @param spOcc data.frame. Object with the coordinate data.
#' @param speciesName character. Name of the species.
#' @param env raster or rasterStack. Environmental data (e.g. typically climatic).
#' @param x name of the field with the coordinate x. Default 'x'
#' @param y name of the field with the coordinate y. Default 'y'
#' @param date name of the field with the values of the reported timestamp of the record. Default NULL
#' @param isoCountry name of the field with the values of the reported country of the record. Default NULL
#' @param classification character. Indicates the thresholds philosophy applied to classify errors in occurrence data. Possible values 'strict','relaxed','custom'
#' @param filterCols logical. Should only the initial input columns be retained in the output (the filtered dataframe)?
#' @return a list of two. First element is a data.frame with profiled occurrence records with their associated profiled labels. Second element is a dataframe with all outputs of the analysis implemented.
#' @note The majority of function parameters can be adjusted but we  provide default values. \cr 
#' We recommend those default values if the user is to use the geospatial data included in the package.\cr
#' but this automatic implementation (occTest + occFilter) missses some analysis to increase speed.
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples \dontrun{
#' #' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#'
#' #load occurrence data
#' occData <- read.csv (system.file('ext/exampleOccData.csv',package = 'occTest'))
#' #load environmental raster
#' renv <- readRDS (system.file('ext/env.rds',package = 'occTest'))
#' #do a simple occTest + occFilter altogether with simplified params
#' outS = occSimpFilter (spOcc = occData,env = renv,speciesName = 'MyFake species', 
#' x = 'MAPX',y = 'MAPY',date = 'DATE',isoCountry = 'COUNTRYRECORD',
#' classification = 'relaxed',
#' filterCols = T)
#' 
#' }
#' @export

#wallace simple wrapping function
occSimpFilter = function(spOcc,env,speciesName='My species',x='x',y='y',
                         date=NULL,isoCountry=NULL,
                         classification='majority',filterCols=T ){
  
  #set up params
  mySettings = occTest::minimalSettings()
  mySettings$tableSettings$x.field <- x
  mySettings$tableSettings$y.field <- y
  mySettings$tableSettings$t.field <- date
  mySettings$tableSettings$c.field <- isoCountry
  
  #run test functions 
  #output is a list of 2 data.frames (the full and  the short, we continue with the full, see below)
  output =  occTest(sp.name = speciesName,
                            sp.table = spOcc, 
                            r.env = env,
                            tableSettings =mySettings$tableSettings)
  
  #select records
  output2 = occTest::occFilter(df = output,errorAcceptance = classification)

  #to consider, maybe you do not want to drag all the columns for Wallace or any other thing, right
  if (filterCols == T) { output2$fitleredDataset[,names (spOcc)]} else {
    output2$fitleredDataset
  }
  
  



  
}
