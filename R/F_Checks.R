#'Checks on the projection of the spatial data
#'
#'Verify that all data are in the same projection
#' @param list.geospatial.objects A list of geospatial objects.Default list includes: 'countries.shapefile','r.env','r.dem','ras.hii','points.proj4string'
#' @return None. Used to generate warning messages.
#' @family checks
#' @examples \dontrun{
#'
#' }

.check.geospatial.data <- function (list.geospatial.objects, verbose=F)  {

                                   
  if (missing (list.geospatial.objects)) {stop ('missing list.geospatial.objects')}
                                    
  if (!is.list (list.geospatial.objects)) {stop ('list.geospatial.objects should be a list')}

  #pf <- parent.frame()
  projections.in <- unlist (lapply (list.geospatial.objects, function(p) {projection(p); return(p)} ))

  if (any(is.na(projections.in))) {
    # id.obj <- which (is.na(projections.in))
    # no.proj.objects <- unlist (list.geospatial.objects[id.obj])
    # stop (print (paste (c('Provide spatial reference data for object:', no.proj.objects),collapse = '  ')))
    stop (print ('One or more geospatial layers do not have projections'))
  }


  tp <- table (unlist (projections.in))
  all.projections.thesame <- as.numeric ((tp)[1] ==length(list.geospatial.objects))
  if (!all.projections.thesame) {
    most.common.proj <- names (tp) [which (tp==max(tp))]
    print ("The most common spatial projection is:")
    print (most.common.proj)

    stop('Please harmonize the coordinate reference system for spatial objects')

  }


}


#checking main fields (inspired in addmainfields from biogeo). I put the number two because it consitutes a version 2 of the functions in biogeo

#'Checking main fields
#'
#'Verify that all main data fields are correctly populated.
#' @param dat A dataframe containing occurrence data for checking.
#' @param xf character. Name of the field where the x coordinate is stored (typically longitude). Default is x.field
#' @param yf character. Name of the field where the y coordinate is stored (typically latitude). Default is y.field
#' @param ef character. Name of the field where the elevation of data collection is stored in the original dataset. Default is e.field.
#' @param tf character. Name of the field where the date of data collection is stored in the original dataset. Default is t.field.
#' @param lf character. Name of the field where the toponim/location of data collection is stored in the original dataset. Default is l.field.
#' @param cf character. Name of the field where the registred country of data collection is stored in the original dataset. Default is c.field.
#' @param idf character. Name of the field of the id of the observation
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @examples \dontrun{
#' @export
#' checkfields()
#' }
.checkfields <- function (dat,
                         xf=x.field,
                         yf=y.field,
                         ef=e.field,
                         tf=t.field,
                         lf=l.field,
                         cf=c.field, 
                         idf = taxonobservation.id,
                         verbose=F){

  nd <- names (dat)
  if ( any (! c(xf,yf) %in% nd )) {stop("need to provide x and y coordinates in your data")}

  #check nulls
  if (is.null (ef) )         {dat$elev <- NA ; if(verbose)  print("Elevation field was NULL. Setting to NA") ; ef<- 'elev' }
  if (is.null (lf) )         {dat$locality <- NA ; if(verbose) print("Locality field was NULL. Setting to NA") ; lf<- 'locality'}
  if (is.null (cf) )         {dat$countryRecorded <- NA ; if(verbose) print("CountryRecorded field was NULL. Setting to NA"); cf<- 'countryRecorded' }
  if (is.null (tf) )         {dat$time <- NA ; if(verbose) print("Time field was NULL. Setting to NA"); tf<- 'time'}
  if (is.null (idf))         {dat$taxonobservationID <- 1:nrow(dat); if (verbose) {print ('taxonObs was NULL. Automatic observation ID implemented')}; idf = 'taxonobservationID'}
  nd <- names (dat)

  #check misspecification
  if (! ef %in% nd )  {stop("elevation field specified not in the occurrence dataframe provided")}
  if (! cf %in% nd )  {stop("contryRecorded field specified not in the occurrence dataframe provided")}
  if (! lf %in% nd )  {stop("locality field specified not in the occurrence dataframe provided")}
  if (! tf %in% nd )  {stop("time field specified not in the occurrence dataframe provided")}
  if (! idf %in% nd )  {stop("ObservatinID field specified not in the occurrence dataframe provided")}
  
  return (dat)

}


#'checking main fields
#'
#'Verify that all main data fields are correctly populated.
#' @param dat A dataframe containing occurrence data for checking.
#' @param xf character. Name of the field where the x coordinate is stored (typically longitude). Default is x.field
#' @param yf character. Name of the field where the y coordinate is stored (typically latitude). Default is y.field
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @note Inspired in addmainfields from biogeo
#' @examples \dontrun{
#' Example<-"goes here"
#' }
#'
.checkdatastr2  <- function (dat,xf=x.field,yf=y.field, verbose=F) {
  cn <- names(dat)
  fn <- c("ID", xf, yf, "Species", "x_original", "y_original",
          "Correction", "Modified", "Exclude", "Reason")
  m <- rep(0, length(fn))
  for (i in 1:length(fn)) {
    m[i] <- match(fn[i], cn, nomatch = 0)
  }
  present <- m > 0
  n <- data.frame(Field = fn, Present = present)
  return(n)
}


#'checking main fields
#'
#'Verify that all main data fields are correctly populated.
#' @param dat A dataframe containing occurrence data for checking.
#' @param species character. Name of the species
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @examples \dontrun{
#' Example<-"goes here"
#' }
#'
.addmainfields2 <- function (dat, species, verbose=F) {

  if (is.na(species)) {
    species <- "Species"
    Species <- NA
  }
  else {
    Species <- dat[, species]
  }
  reqNames <- c("ID", species, "x_original", "y_original",
                "Correction", "Modified", "Exclude", "Reason")
  missingNames <- reqNames[!sapply(reqNames, FUN = function(x) x %in%
                                     names(dat))]
  z <- data.frame(dat, ID = 1:nrow(dat), Species, x = NA, y = NA,
                  x_original = NA, y_original = NA, Correction = "........",
                  Modified = Sys.time(), Exclude = 0, Reason = "........",
                  stringsAsFactors = F)
  z <- z[, c(names(dat), missingNames)]
  return(z)
}



#'checking main fields
#'
#'Check for missing fields.
#' @param dat A dataframe containing occurrence data for checking.
#' @param fields character. Fields to check
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @examples \dontrun{
#' Example<-"goes here"
#' }
# fieldsMissing2 <- function (dat, fields) {
#   ck <- checkDataStr2(dat,xf = x.field, yf = y.field)
#   if (any(ck$Present == FALSE) == TRUE) {
#     fa <- ck[ck$Present == FALSE, 1]
#     fd <- as.character(fa)
#     msg <- paste(fd, collapse = ", ")
#     stop(paste("Fields missing: ", msg))
#   }
# }


#'STATUS TRACKER AND WRITE OUTPUT
#'
#'Track status and write useful output
#' @param dataset.to.continue A dataframe containing occurrence data for checking.
#' @param wfo write full output
#' @param wso write simple output
#' @param xf The dataframe field containing the x values (e.g. "longitude")
#' @param yf The dataframe field containing the y values (e.g. "latitude")
#' @param od The output directory to use
#' @param obf Output base filename
#' @param sp character. Name of the species
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @examples \dontrun{
#' Example<-"goes here"
#' }
.status.tracker.and.escaping <- function (dataset.to.continue=dat,
                                         wfo=write.full.output,
                                         wso=write.simple.output,
                                         xf=x.field,
                                         yf=y.field,
                                         od=output.dir,
                                         obf=output.base.filename,
                                         sp=sp.name, verbose=F){

  if (nrow (dataset.to.continue) != 0) {return(NULL)}
  if (nrow (dataset.to.continue) == 0) {print ('Workflow finished')}
  if (nrow (dataset.to.continue) == 0) {

    ### exit control flow
    pf <- parent.frame()
    all.potential.output.qdf  <- grep(pattern = 'dat.Q',ls(name = pf),value = T)

    dat.out.list <- lapply (all.potential.output.qdf,
                            function (i){if (exists(i,pf)){get (i,pf)}})

    if (length(dat.out.list) == 1) {full.qaqc <- unlist(dat.out.list)}
    if (length(dat.out.list) > 1)  {full.qaqc <- Reduce (plyr::rbind.fill, dat.out.list)}

    full.qaqc$qualifiers<- NA
    full.qaqc$quality.label<- full.qaqc$quality.grade

    short.qaqc <-full.qaqc[, c('ID',xf,yf,'quality.grade','qualifiers','quality.label')]


    if (wfo==T) {
      write.csv (full.qaqc,   paste0(od,'/',obf,'_',sp,'_long.csv'),row.names = F)
    }

    if (wso==T) {
      write.csv (short.qaqc,   paste0(od,'/',obf,'_',sp,'_short.csv'),row.names = F)
    }

    output.function <- list (occ_full_profile=full.qaqc, occ_short_profile=short.qaqc)
    return (output.function)
  }

}


