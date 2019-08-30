#### FUNCTIONS FOR INVASIVE AND NATIVE RANGES AND COUNTRY STATUS

#' @title Check native and alien ranges in countries
#'
#' @description performs several queries to 
#' @details right now based on the package originatr. Future implementation may query more databases.
#'
#' @param spName character. Species name in the form of Genus species
#' @param resolveNative logical. Should the function attempt to find countries where species are considered native?
#' @param resolveAlien logical. Should the function attempt to find countries where species are considered introduced?
#' @param verbose logical. Want to print information during the process?
#' @return list with two vectors: ntvCtry and invCtry, showing the countries in ISO3 standard codes
#' @keywords internal
#'
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
#' 

# 
# nsr("Pinus ponderosa", country = "United States")
# is_native
# library (countrycode)
# 
# countrycode('spain', 'country.name', 'iso3c')
# 
# devtools::install_github("ropensci/originr")
# 
# 
# 
# ### examples start here
# 
# spName = 'Ailanthus altissima'
# 
# if (interactiveMode & is.null (ntv.range) ){
#   resolveNative <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species native range. Do you want to infer from global databases?")
# }
# if (interactiveMode & is.null (inv.range) ){
#   resolveAlien <- if (interactive())  askYesNo(default = F,msg = "You have not provided countries for the species alien range. Do you want to infer from global databases?")
# }
# if (resolveAlienNative) {resolveAlien = T ;resolveNative = T }
# 
# if (any (resolveAlien, resolveNative)){
#   



nativeStatusCtry <- function (spName,resolveNative=T,resolveAlien=T, verbose=T){

  #check data from flora Europaea
  try (floraEurope <- originr::flora_europaea(spName), silent=T)
  
  #check gisd (it also contains info on native range)
  try (gisdInvasive <- originr::gisd(x = spName), silent=T)
  
  #check info en several invasive checklists
  if (resolveAlien){
    #not implemented
    #try (eolInvasive <- originr::eol(name = spName),silent=T)
    try (griisInvasive <-originr::griis(name = spName), silent = T)
  }
  
  #check data from BIEN native species resolver 
  data (nsr_countries,package = 'originr')
  listNSR = lapply (nsr_countries, function (x,s=spName){
    try (originr::nsr(s,x),silent = T)
  })
  try (listNsrClean <- lapply (listNSR, function (x) if (class(x)=='try-error' ) {NULL} else {x}), silent=T )
  try (listNsrClean <- do.call (rbind,listNsrClean), silent= T)
  
  #gather data for natives
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
  
  #gather data for alien ranges in countries
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
  ntvCtryResolved <- occProfileR:::ctryToIso3(ntvCtryResolved)
  invCtryResolved <- occProfileR:::ctryToIso3(invCtryResolved)
  if (any(ntvCtryResolved %in% invCtryResolved)) {
    # print (paste ('Species is considered native and invasive in',paste(ntvCtryResolved[ntvCtryResolved %in% invCtryResolved],collapse = ','),sep = ': ' ) )
    # print ('We will consider it only as in the native range')
    # 
    invCtryResolved <- invCtryResolved[! invCtryResolved %in% ntvCtryResolved]
    
  }
  out <- list (ntvCtry= unique (ntvCtryResolved), invCtry= unique(invCtryResolved))
  return (out)
  
}






#' @title Convert country names to ISO3 codes
#'
#' @description using package countrycode 
#' @details right now not implemented with fuzzy matching, but is case insensitive.
#'
#' @param x character. country name
#' @keywords internal
#'
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }




ctryToIso3 <- function (x){ 

  if (!is.character(x) & !is.null(x)) {stop("Country should be character")}
  if (length(x)==0 | is.null(x) ) {outIso3<- NULL}
  if (length(x)>0 ) {
    #try to correCt for potential bad spelling
    template<-GNRS::GNRS_template(nrow = length(x))
    template$country<- x
    gnrsResults  = GNRS:::GNRS(political_division_dataframe = template)
    iso2Ctry <- unlist(gnrsResults$country_iso)
    
    #get iso3 codes
    outIso3 <- as.character (sapply (iso2Ctry, function (x){countrycode::countrycode(x, 'iso2c', 'iso3c')},simplify=T))
    outIso3 <- outIso3[!is.na(outIso3)]
    if(length(outIso3) == 0) {outIso3 <- NULL}
    
  } 
  return (outIso3)
  
  
  
  
}
  

