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
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family spStatus
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export
nativeStatusCtry <- function (spName,xydat, resolveNative=T,resolveAlien=T, verbose=T){


  #CHECK FLORA EUROPAEA
  try (floraEurope <- originr::flora_europaea(spName,messages = verbose), silent=T)
  
  #CHECK GISD (it also contains info on native range)
  try (gisdInvasive <- originr::gisd(x = spName,messages = verbose), silent=T)
  
  #CHECK BIEN NRS (contains both invasive and native )
  ctryFromCoords =  occTest::.coords2country (xydat)
  ctryFromCoords =  unique (ctryFromCoords)
  ctry4GNRS = countrycode::countrycode(sourcevar = ctryFromCoords,origin = 'iso3c',destination = 'country.name')
  #we apply over cuz GNRS only accepts one query at at ime
  listNSR = lapply (ctry4GNRS, function (x,s=spName){
    try (originr::nsr(s,x),silent = T)
  })
  listNsrClean = try (expr = Filter(Negate(is.null), listNSR), silent=T)
  listNsrClean = try (expr = do.call (rbind,listNsrClean)    , silent=T)
  
  #CHECK INFO INVASIVE LISTS 
  if (resolveAlien){
    #not implemented
    #try (eolInvasive <- originr::eol(name = spName),silent=T)
    try (griisInvasive <-originr::griis(name = spName), silent = T)
  }
  
  #GATHER DATA NATIVES
  if (resolveNative){
    ntvCtryResolved <- NULL
    if (exists('floraEurope') & class(floraEurope)=="list") {ntvCtryResolved <- c(ntvCtryResolved,floraEurope$native)}
    if (exists('listNsrClean') & class(listNsrClean)=="data.frame") {
      ntvNSR <- listNsrClean [listNsrClean$native_status == 'N', 'country']
      if (length(ntvNSR)==0) {ntvNSR<- NULL}
      ntvCtryResolved <- c(ntvCtryResolved, ntvNSR)
    }
    if (exists('gisdInvasive') & class(gisdInvasive)=="list") {
      ntvCtryResolved <- c (ntvCtryResolved,gisdInvasive[[1]]$native_range)
    }
    
    ntvCtryResolved <- ntvCtryResolved[!is.na(ntvCtryResolved)]
    if (length(ntvCtryResolved)==0 & verbose) {"no native countries found. Native countries set to NULL"}
    if (length(ntvCtryResolved)>0 & verbose) {paste("Native countries found",ntvCtryResolved)}
  }
  
  #GATHER DATA INVASIVE 
  if (resolveAlien){
    invCtryResolved <- NULL
    if (exists('floraEurope') & class(floraEurope)=="list") {invCtryResolved <- c(invCtryResolved,floraEurope$exotic)}
    if (exists('gisdInvasive') & class(gisdInvasive)=="list") {
      invCtryResolved <- c (invCtryResolved,gisdInvasive[[1]]$alien_range)
    }
    if (exists('griisInvasive') & class(griisInvasive)=="data.frame") {
      griisInvasive <- griisInvasive[griisInvasive$Origin=='Alien',]
      if (nrow(griisInvasive)>0){
        invCtryResolved <- c (invCtryResolved,griisInvasive$Country)
      }
      
    }
    if (exists('listNsrClean') & class(listNsrClean)=="data.frame") {
      invNSR <- listNsrClean [listNsrClean$native_status %in% c('I','A'), 'country']
      if (length(invNSR)==0) {invNSR<- NULL}
      invCtryResolved <- c(invCtryResolved, invNSR)
    }
    
    
    
    invCtryResolved <- invCtryResolved[!is.na(invCtryResolved)]
    if (length(invCtryResolved)==0 & verbose) {"no Alien ranges found in countries. Native countries set to NULL"}
    if (length(invCtryResolved)>0 & verbose) {paste("Alien ranges found in countries",paste0 (invCtryResolved,collapse =';'))}
    
  }


  #convert country names to iso3 country codes
  ntvCtryResolved <- occTest::ctryToIso3(ntvCtryResolved)
  invCtryResolved <- occTest::ctryToIso3(invCtryResolved)
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
#' @details right now not implemented with fuzzy matching, but is case insensitive. Tow methods implemented, 'countrycode' and 'GNRS'
#' @param x character. country name
#' @param method character. Package name used to derive IS03 codes. Options are 'countrycode' (default) or 'GNRS'.
#' @keywords spStatus
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family spStatus
#' @examples \dontrun{
#' example<-"goes here"
#' }
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
    outIso3 <- as.character (sapply (iso2Ctry, function (x){countrycode::countrycode(x, 'iso2c', 'iso3c')},simplify=T))
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
