# .check.geospatial.dat ====
#' @title Checks on the projection of the spatial data
#' @description Verify that all data are in the same projection
#' @param list.geospatial.objects A list of geospatial objects.Default list includes: 'countries.shapefile','r.env','r.dem','ras.hii','points.proj4string'
#' @param verbose logical. Print messages?
#' @author Josep M Serra Diaz
#' @return None. Used to generate warning messages.
#' @family checks
#' @examples {
#' example<-"goes here"
#' }

.check.geospatial.data <- function (list.geospatial.objects, verbose=FALSE)  {

                                   
  if (missing (list.geospatial.objects)) {stop ('missing list.geospatial.objects')}
                                    
  if (!is.list (list.geospatial.objects)) {stop ('list.geospatial.objects should be a list')}

  #pf <- parent.frame()
  projections.in <- unlist (lapply (list.geospatial.objects, function(p) {raster::projection(p); return(p)} ))

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


# .checkfields ====
#' @title Checking main fields
#' @descriptoin Verify that all main data fields are correctly populated.
#' @details  checking main fields (inspired by \link[biogeo]{addmainfields} . 
#' @param dat A dataframe containing occurrence data for checking.
#' @param xf character. Name of the field where the x coordinate is stored (typically longitude). Default is x.field
#' @param yf character. Name of the field where the y coordinate is stored (typically latitude). Default is y.field
#' @param ef character. Name of the field where the elevation of data collection is stored in the original dataset. Default is e.field.
#' @param tf character. Name of the field where the date of data collection is stored in the original dataset. Default is t.field.
#' @param lf character. Name of the field where the toponim/location of data collection is stored in the original dataset. Default is l.field.
#' @param cf character. Name of the field where the registered country of data collection is stored in the original dataset. Default is c.field.
#' @param idf character. Name of the field of the id of the observation
#' @param verbose logical. Print messages? Default to FALSE
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @author Mark Robertson and Vernon Visser (original function), Josep M Serra Diaz (modifs)
.checkfields <- function (dat,
                         xf,
                         yf,
                         ef,
                         tf,
                         lf,
                         cf, 
                         idf,
                         verbose=FALSE){

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
  if (! cf %in% nd )  {stop("countryRecorded field specified not in the occurrence dataframe provided")}
  if (! lf %in% nd )  {stop("locality field specified not in the occurrence dataframe provided")}
  if (! tf %in% nd )  {stop("time field specified not in the occurrence dataframe provided")}
  if (! idf %in% nd )  {stop("ObservationID field specified not in the occurrence dataframe provided")}
  
  return (dat)

}


# .checkdatastr2 ====
#' @title Check data structure
#' @description Verify that all main data fields are correctly structured
#' @details Inspired by \link[biogeo]{checkdatastr} but modified (hence number 2 after the function original name)
#' @param dat A dataframe containing occurrence data for checking.
#' @param xf character. Name of the field where the x coordinate is stored (typically longitude). Default is x.field
#' @param yf character. Name of the field where the y coordinate is stored (typically latitude). Default is y.field
#' @param verbose logical. Print messages? Defaults to FALSE
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @author Mark Robertson and Vernon Visser (original function), Josep M Serra Diaz (modifs)

.checkdatastr2  <- function (dat,xf,yf, verbose=FALSE) {
  cn <- names(dat)
  fn <- c("roworder", xf, yf, "Species", "x_original", "y_original",
          "Correction", "Modified", "Exclude", "Reason")
  m <- rep(0, length(fn))
  for (i in 1:length(fn)) {
    m[i] <- match(fn[i], cn, nomatch = 0)
  }
  present <- m > 0
  n <- data.frame(Field = fn, Present = present)
  return(n)
}


# .addmainfields2 ====

#' @title  Add main fields 
#' @description Incorporate fields in the initial data frame
#' @details Inspired by \link[biogeo]{addmainfields} but modified (hence number 2 after the function original name)
#' @param dat A dataframe containing occurrence data for checking.
#' @param species character. Name of the species
#' @param verbose logical. Print messages? Default FALSE
#' @return Original dataframe, dat. 
#' @author Mark Robertson and Vernon Visser (original function), Josep M Serra Diaz (modifs)
#' @family checks

.addmainfields2 <- function (dat, species, verbose=FALSE) {

  if (is.na(species)) {
    species <- "Species"
    Species <- NA
  }
  else {
    Species <- dat[, species]
  }
  reqNames <- c("roworder", species, "x_original", "y_original",
                "Correction", "Modified", "Exclude", "Reason")
  equal <- reqNames[sapply(reqNames, FUN = function(x) x %in%
                                     names(dat))]
  
  if (length(equal)>1) {warning(paste("Table fields: ",equal,"already existed in your table, not overwritten by fieldchecks.",'\n',"check consistency of the field meanings in your table with occTest"))}
  missingNames <- reqNames[!sapply(reqNames, FUN = function(x) x %in%
                                     names(dat))]
  newdf = lapply(missingNames, function(n) {data.frame (a=rep(NA,length.out=nrow(dat)))})
  newdf = do.call (what = cbind, newdf)
  names(newdf) = missingNames
  if ('Exclude' %in% missingNames) newdf$Exclude = 0
  if ('roworder' %in% missingNames) newdf$roworder = 1:nrow (newdf)
  
  z <- data.frame(dat, newdf,stringsAsFactors = FALSE)
  return(z)
}

# .status.tracker.and.escaping  ====
#' @title Workflow status tracker
#' @description Track status and write useful output
#' @param dataset.to.continue A dataframe containing occurrence data for checking.
#' @param wfo write full output
#' @param wso write simple output
#' @param xf The dataframe field containing the x values (e.g. "longitude")
#' @param yf The dataframe field containing the y values (e.g. "latitude")
#' @param rsd logical. Return spatial data?
#' @param as list. analysis settings
#' @param ws list. writoutsettings
#' @param ts tlist. able setting
#' @param od The output directory to use
#' @param obf Output base filename
#' @param sp character. Name of the species
#' @param verbose logical. Print messages? Defaults to FALSE
#' @return Original dataframe, dat.  Used primarily to generate warning messages.
#' @family checks
#' @author Josep M Serra Diaz

.status.tracker.and.escaping <- function (dataset.to.continue,
                                         wfo,
                                         wso,
                                         xf,
                                         yf,
                                         od,
                                         rsd ,
                                         obf,
                                         sp, verbose=FALSE, 
                                         as,
                                         ws,ts){

  if (nrow (dataset.to.continue) != 0) {return(NULL)}
  if (nrow (dataset.to.continue) == 0) {print ('Workflow finished')}
  if (nrow (dataset.to.continue) == 0) {

    ### exit control flow
    pf <- parent.frame()
    all.potential.output.qdf  <- grep(pattern = 'dat.Q',ls(name = pf),value = TRUE)

    dat.out.list <- lapply (all.potential.output.qdf,
                            function (i){if (exists(i,pf)){get (i,pf)}})

    if (length(dat.out.list) == 1) {full.qaqc <- unlist(dat.out.list)}
    if (length(dat.out.list) > 1)  {full.qaqc <- Reduce (plyr::rbind.fill, dat.out.list)}

    #write outputs
    if(wfo){
      sp2 =  .join.spname(sp)
      newdir = paste0(od,'/',sp2)
      dir.create(newdir,recursive = TRUE,showWarnings = FALSE)
      written = try(utils::write.csv(full.qaqc,  
                              paste0(newdir,'/',obf,
                                     '_',sp,'_long.csv'),
                              row.names = FALSE),silent = TRUE)
      if(inherits(written,'try-error')) save(list = 'full.qaqc',file = paste0(newdir,'/',obf,'_',sp,'_long.RData'))
      if(inherits(written,'try-error')) try(file.remove(paste0(newdir,'/',obf,'_',sp,'_long.csv')), silent=TRUE )
    }
    

    #output.function = list(occTest_full=full.qaqc, occTest_short=short.qaqc)
    output.function = full.qaqc

    attr(output.function,"class")<-c("occTest",class(output.function))
    
    if(!rsd){
      as$countryStatusRange$countries.shapefile<-NULL
      as$humanDetection$ras.hii<-NULL
      as$humanAnalysis$methodHyperHumanDetection<-NULL
      as$rangeAnalysis$countries.shapefile<-NULL
    }
    
    
    attr(output.function,"Settings")<-list(tableSettings=ts,analysisSettings=as,writeoutSettings=ws)
    
    
    return(output.function)
    
    
    

    

  
  }

}


