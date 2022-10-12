### Functions related to identifying native or invasive status ###

# nativeStatusCtry ====
#' @title Check native and alien ranges in countries
#' @description performs several queries based on the orignir package to identify native and invasive range status for a given species
#' @details right now based on the package originr, which queries on Flora Europea, GISD or Native species resolver. Future implementation may query more databases.
#' @param spName character. Species name in the form of Genus species
#' @param xydat dataframe Species longitude latitude coordinates.
#' @param resolveNative logical. Should the function attempt to find countries where species are considered native?
#' @param resolveAlien logical. Should the function attempt to find countries where species are considered introduced?
#' @param verbose logical. Want to print information during the process?
#' @return list with two vectors: ntvCtry and invCtry, showing the countries in ISO3 standard codes
#' @keywords internal
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @family spStatus
#' @notes originr is a non-standard package. \cr
#' Working on updates and alternative workflows as the package has been abandonded \cr
#' source code here: https://github.com/ropensci-archive/originr
#' @examples \donttest{
#' #needs spocc library
#' df <- spocc::occ(query = 'Pseudotsuga menziesii')
#' occ.data <- spocc::occ2df(df)
#' nativeStatusCtry(spName = 'Pseudotsuga menziesii',xydat = occ.data[,2:3])
#' }
#' @export
nativeStatusCtry <- function (spName,xydat, resolveNative=TRUE,resolveAlien=TRUE, verbose=TRUE){


  #CHECK FLORA EUROPAEA
  try (floraEurope <- originr::flora_europaea(spName,messages = verbose), silent=TRUE)
  
  #CHECK GISD (it also contains info on native range)
  try (gisdInvasive <- originr::gisd(x = spName,messages = verbose), silent=TRUE)
  
  #CHECK BIEN NRS (contains both invasive and native )
  ctryFromCoords =   .coords2country (xydat)
  ctryFromCoords =  unique (ctryFromCoords)
  ctry4GNRS = countrycode::countrycode(sourcevar = ctryFromCoords,origin = 'iso3c',destination = 'country.name')
  #we apply over cuz GNRS only accepts one query at at ime
  listNSR = lapply (ctry4GNRS, function (x,s=spName){
    try (originr::nsr(s,x),silent = TRUE)
  })
  listNsrClean = try (expr = Filter(Negate(is.null), listNSR), silent=TRUE)
  listNsrClean = try (expr = do.call (rbind,listNsrClean)    , silent=TRUE)
  
  #CHECK INFO INVASIVE LISTS 
  if (resolveAlien){
    #not implemented
    #try (eolInvasive <- originr::eol(name = spName),silent=TRUE)
    try (griisInvasive <-originr::griis(name = spName), silent = TRUE)
  }
  
  #GATHER DATA NATIVES
  if (resolveNative){
    ntvCtryResolved <- NULL
    if (exists('floraEurope') & inherits(floraEurope,"list")) {ntvCtryResolved <- c(ntvCtryResolved,floraEurope$native)}
    if (exists('listNsrClean') & inherits(listNsrClean,"data.frame")) {
      ntvNSR <- listNsrClean [listNsrClean$native_status == 'N', 'country']
      if (length(ntvNSR)==0) {ntvNSR<- NULL}
      ntvCtryResolved <- c(ntvCtryResolved, ntvNSR)
    }
    if (exists('gisdInvasive') & inherits(gisdInvasive,'list')) {
      ntvCtryResolved <- c (ntvCtryResolved,gisdInvasive[[1]]$native_range)
    }
    
    ntvCtryResolved <- ntvCtryResolved[!is.na(ntvCtryResolved)]
    if (length(ntvCtryResolved)==0 & verbose) {"no native countries found. Native countries set to NULL"}
    if (length(ntvCtryResolved)>0 & verbose) {paste("Native countries found",ntvCtryResolved)}
  }
  
  #GATHER DATA INVASIVE 
  if (resolveAlien){
    invCtryResolved <- NULL
    if (exists('floraEurope')) {if (inherits(floraEurope,"list")) {invCtryResolved <- c(invCtryResolved,floraEurope$exotic)}}  
    if (exists('gisdInvasive')) {if (inherits(gisdInvasive,"list")) {
      invCtryResolved <- c (invCtryResolved,gisdInvasive[[1]]$alien_range)
    }
      }  
    if (exists('griisInvasive'))  {if (inherits(griisInvasive,"data.frame")) {
      griisInvasive <- griisInvasive[griisInvasive$Origin=='Alien',]
      if (nrow(griisInvasive)>0){
        invCtryResolved <- c (invCtryResolved,griisInvasive$Country)
      }
      
    }}
    if (exists('listNsrClean'))  {if (inherits(listNsrClean,"data.frame")) {
      invNSR <- listNsrClean [listNsrClean$native_status %in% c('I','A'), 'country']
      if (length(invNSR)==0) {invNSR<- NULL}
      invCtryResolved <- c(invCtryResolved, invNSR)
    }}
    
    
    
    invCtryResolved <- invCtryResolved[!is.na(invCtryResolved)]
    if (length(invCtryResolved)==0 & verbose) {"no Alien ranges found in countries. Native countries set to NULL"}
    if (length(invCtryResolved)>0 & verbose) {paste("Alien ranges found in countries",paste0 (invCtryResolved,collapse =';'))}
    
  }


  #convert country names to iso3 country codes
  ntvCtryResolved <-  ctryToIso3(ntvCtryResolved)
  invCtryResolved <-  ctryToIso3(invCtryResolved)
  if (any(ntvCtryResolved %in% invCtryResolved)) {
    # print (paste ('Species is considered native and invasive in',paste(ntvCtryResolved[ntvCtryResolved %in% invCtryResolved],collapse = ','),sep = ': ' ) )
    # print ('We will consider it only as in the native range')
    # 
    invCtryResolved <- invCtryResolved[! invCtryResolved %in% ntvCtryResolved]
    
  }
  out <- list (ntvCtry= unique (ntvCtryResolved), invCtry= unique(invCtryResolved))
  return (out)
  
}





# ctryToIso3 ====
#' @title Convert country names to ISO3 codes
#' @description Froma character it uses different methods to derive country ISO3 digit codes
#' @details right now not implemented with fuzzy matching, but is case insensitive. Methods implemented, 'countrycode' and 'GNRS'
#' @param x character. country name
#' @param method character. Package name used to derive IS03 codes. Options are 'countrycode' (default) or 'GNRS'.
#' @keywords spStatus
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr) (adpatation) \cr 
#' \link[countrycode]{countrycode-package} Vincent Arel-Bundock vincent.arel-bundock@@umontreal.ca ,\cr
#' \link[GNRS]{GNRS}  Brad Boyle, Brian Maitner 
#' @family spStatus
#' @seealso \link[countrycode]{countrycode} and \link[GNRS]{GNRS}

ctryToIso3 <- function (x,method='countrycode'){ 
  
  #initial checks
  if (!is.character(x) & !is.null(x)) {stop("Country should be character")}
  if (length(x)==0 | is.null(x) ) {return (NULL)}
  
  
  if (method=='GNRS' ) {
    #try to correCt for potential bad spelling
    tableTemplate<-GNRS::GNRS_template(nrow = length(x))
    tableTemplate$country<- x
    gnrsResults  = GNRS::GNRS(political_division_dataframe = tableTemplate)
    iso2Ctry <- unlist(gnrsResults$country_iso)
    
    #get iso3 codes
    outIso3 <- as.character (sapply (iso2Ctry, function (x){countrycode::countrycode(x, 'iso2c', 'iso3c')},simplify=TRUE))
    outIso3 <- outIso3[!is.na(outIso3)]
    if(length(outIso3) == 0) {outIso3 <- NULL}
    
  } 
  
  if (method=='countrycode' ) {
    
    #direct matching
    outIso3 = countrycode::countrycode (sourcevar = x, origin = 'country.name',destination = 'iso3c')
    
    #fuzzy matching  using stringdist package (not implemented yet)
    # codesAll = countrycode::codelist
    # colid  = grep (pattern = 'country.name.en.regex',names(codesAll))
    # colid  = grep (pattern = 'cldr.variant.es',names(codesAll))
    # dfForMatch = codesAll[,colid]
    # adist (x,dfForMatch,method='osa',maxDist = 50)
    # 
    # apply (dfForMatch,MARGIN = 1,FUN = function (n){adist(n,x[1])})
    #
    
    
    
  } 
  
  return (outIso3)
  

}
