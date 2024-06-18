# centroid_assessment_BIEN ====
#' @title Identify occurrences via BIEN proximity analysis
#' @author Brian Maitner & Pep Serra-Diaz
#' @description
#' This function is designed to use the BIEN centroid detection method which provides both relative and absolute distances to know centroids.  The lower the relative distance, the more likely that the occurrence is a centroid.
#' @param occurrences A occurrence sf data frame containing the correct fields. See example for more details.
#' @param centroid_data Centroid sf data.frame from BIEN database
#' @param relative_distance_threshold Threshold below which coordinates are assumed to be centroids. Default value is...
#' @import dplyr
#' @importFrom tibble remove_rownames
#' @note
#' Relative centroid distance is calculated as: (distance between the occurrence record and centroid)/(distance between centroid and the furthest from it within the region)
#' @return A data.frame containing information about centroid inference.  
#' is_centroid: Logical.  1 = Yes, 0 = No, NA = assessment could not be made.
#' centroid_type: The type of centroid. Main refers to centroid of the largest polygon.  BB refers to the polygon bounding box.
#' relative_distance: Relative distance between the occurrence point and the nearest centroid (see note)
#' centroid_distance: Distance (decimal degrees) between the occurrence point and the nearest centroid.
#' dist_max: Furthest distance from the centroid within the political division. 
#' @export
#' @examples
#' \dontrun{
#' # Load required libraries
#'
#' library(occTest)
#' library(piggyback)
#' library(sf)
#' library(tidyverse)
#'
#'# Create temporary directory to store the centroid data
#'
#'temp_dir <- tempdir()
#'
#'# Download the data from Github
#'
#'pb_download(file = "centroids_long_format.gpkg",
#'            repo = "pepbioalerts/occTest",
#'            tag = "data",
#'            dest = temp_dir,
#'            overwrite = TRUE)
#'
#'# Load centroid file
#'
#'centroids <- st_read(file.path(temp_dir,"centroids_long_format.gpkg"))
#'
#'# Load example data set
#'
#'occData <- read.csv(system.file('ext/exampleOccData.csv',
#'                                package = 'occTest')) |>
#'  filter(!is.na(MAPX) & !is.na(MAPY) ) |>
#'  st_as_sf(coords = c("MAPX", "MAPY")) |>
#'  mutate(country  = COUNTRYRECORD)
#'
#'
#'# Run function and append to existing data
#'
#'occData |>
#'  bind_rows(  centroid_assessment_BIEN(occurrences = occData,
#'                                       centroid_data = centroids))
#'}                                       
centroid_assessment_BIEN <- function(occurrences,
                                     centroid_data,
                                     relative_distance_threshold = 0.001){

  #checks for data format

    if(!inherits(x = occurrences, what = "sf")){
      stop("Occurrences should be supplied as an sf data.frame")
    }
    
    if(!inherits(x = centroid_data, what = "sf")){
      stop("Occurrences should be supplied as an sf data.frame")
    }
  
  #add unique ID field

    occurrences <-
      occurrences |>
      dplyr::mutate(cent_ID = 1:dplyr::n())
  
  ### NEW
  #get info country from occurrences (whether inferred or original country record)
    occurrences$countryBIEN <- ifelse (is.na (occurrences$country),occurrences$countryCoordinates,occurrences$country)
    occurrences$countryBIEN_commment <- ifelse (is.na (occurrences$country),'region inferred from coordinates','region inferred from country recorded')

  #Add empty country and state columns if needed  
    
    if(!"state_province" %in% colnames(occurrences)){
      occurrences$state_province <- NA
    }
    
    if(!"county_parish" %in% colnames(occurrences)){
      occurrences$county <- NA
    }
    
  # check that country is an ISO3 code
    
    if(all(!((nchar(occurrences$countryBIEN) == 3) | (is.na(occurrences$countryBIEN))))){
      stop("Fields country or countryBIEN should be ISO3 codes")
    }
  
    
  #country data
    country_centroids <- NULL
  
  ### NEW :see that I changed the field to countryBIEN 
  #(not changed all the way through tho, double check)
  for(i in 1:length(unique(occurrences$countryBIEN))){
    
    # If there are no countries, move on
    if(length(unique(occurrences$countryBIEN)) < 1){break}
    
    country_i <- unique(occurrences$countryBIEN)[i]
    
    # If the country is NA, move on
    if (is.na(country_i) | country_i == ""){ next()}
    
    centroid_i <-
      centroid_data |>
      dplyr::filter(gid_0 == country_i &
                is.na(gid_1) &
                is.na(gid_2) )
    
    # If not centroid information, move on
    if (nrow (centroid_i)==0){ next()}
    
    occurrences_i <-
      occurrences |>
      dplyr::filter(countryBIEN == country_i &
                is.na(state_province) &
                is.na(county) )
    
    sf::st_crs(occurrences_i) <- st_crs(NULL)
    sf::st_crs(centroid_i) <- st_crs(NULL)
    
    if(nrow(centroid_i) > 0){
      
      #set the crs to null so that the units will be decimal degrees (to match the max distances used in BIEN)
      
        sf::st_crs(occurrences_i) <- sf::st_crs(NULL)
        sf::st_crs(centroid_i) <- sf::st_crs(NULL)
      
      #calculate distances  
      
        dists_dd <- sf::st_distance(occurrences_i, centroid_i)
      
      #calculate minimum distances and get the closest centroid
        
        dists_dd |>
          apply(MARGIN = 1, FUN = which.min) -> dist_dd_index
        
        dists_dd |>
          apply(MARGIN = 1, FUN = min) -> mins_dd
      
        closest_centroids_i <- centroid_i[dist_dd_index,]
      
      closest_centroids_i <-
        closest_centroids_i |>
        dplyr::mutate(cent_ID = occurrences_i$cent_ID,
                      centroid_distance = mins_dd,
               relative_distance = centroid_distance/dist_max) |>
        sf::st_drop_geometry()
      
      if(nrow(closest_centroids_i) != nrow(occurrences_i)){
        stop("occurrence vs centroid mismatch in countries") }
      
      country_centroids |>
        dplyr::bind_rows(closest_centroids_i) -> country_centroids
      
      
    }else{
      
      #stop("Do something if there aren't any centroids (NAs prob)")
      
      occurrences_i |>
        dplyr::select(.data,cent_ID,country,state_province,county)|>
        dplyr::rename(county_parish=county)|>
        sf::st_drop_geometry(.data)|>
        dplyr::bind_rows(country_centroids) -> country_centroids
      
      
    }
    
    
  }#i country loop
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #state data
  
  state_centroids <- NULL
  
  states <-
    occurrences |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(country) &
             !is.na(state_province)&
             is.na(county))|>
    dplyr::select(country,state_province)|>
    unique()
  
  
  for(i in 1:nrow(states)){
    
    if(nrow(states)<1){break}
    
        
    country_i <-states[i,'country']
    state_i <-states[i,'state_province']
    
    centroid_i <-
      centroid_data |>
      dplyr::filter( country == country_i &
                       state_province==state_i &
                is.na(county_parish) )
    
    occurrences_i <-
      occurrences |>
      dplyr::filter( country == country_i &
                       state_province==state_i &
                is.na(county) )
    
    sf::st_crs(occurrences_i) <- st_crs(NULL)
    sf::st_crs(centroid_i) <- st_crs(NULL)
    
    if(nrow(centroid_i) > 0){
      
      # dists_m <- st_distance(occurrences_i, centroid_i)
      # 
      # dists |>
      #   apply(MARGIN = 1, FUN = which.min) -> dist_m_index
      # 
      # dists |>
      #   apply(MARGIN = 1, FUN = min) -> mins_m
      
      
      #set the crs to null so that the units will be decimal degrees (to match the max distances used in BIEN)
      
      sf::st_crs(occurrences_i) <- sf::st_crs(NULL)
      sf::st_crs(centroid_i) <- sf::st_crs(NULL)
      
      #calculate distances  
      
      dists_dd <- sf::st_distance(occurrences_i, centroid_i)
      
      dists_dd |>
        apply(MARGIN = 1, FUN = which.min) -> dist_dd_index
      
      dists_dd |>
        apply(MARGIN = 1, FUN = min) -> mins_dd
      
      closest_centroids_i <- centroid_i[dist_dd_index,]
      
      closest_centroids_i <-
        closest_centroids_i |>
        dplyr::mutate(cent_ID = occurrences_i$cent_ID,
               centroid_distance = mins_dd,
               relative_distance = centroid_distance/dist_max)|>
        sf::st_drop_geometry()
      
      if(nrow(closest_centroids_i) != nrow(occurrences_i)){
        stop("occurrence vs centroid mismatch in states")}
      
      state_centroids |>
        dplyr::bind_rows(closest_centroids_i) -> state_centroids
      
      
    }else{
      
      occurrences_i |>
        dplyr::select(cent_ID,country,state_province,county)|>
        dplyr::rename(.data,county_parish=county)|>
        sf::st_drop_geometry(.data)|>
        dplyr::bind_rows(state_centroids) -> state_centroids
      
      
    }
    
    
    
  }#i state loop
  
  #check lengths
  occurrences |>
    dplyr::filter(!is.na(country) &
             !is.na(state_province)&
             is.na(county))|>
    sf::st_drop_geometry()->state_occs
  
  if(nrow(state_occs) > 0 ){
    
    if(nrow(state_occs) != nrow(state_centroids)){stop("mismatch in state centroid lengths")}  
    
  }
  
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # county data
  
  county_centroids <- NULL
  
  counties <-
    occurrences |>
    sf::st_drop_geometry(.data) |>
    dplyr::filter(!is.na(country) &
             !is.na(state_province)&
             !is.na(county))|>
    dplyr::select(country,state_province,county)|>
    unique()
  
  for(i in 1:nrow(counties)){
    
    if(nrow(counties) <1 ){ break }
    
    country_i <- counties[i,'country']
    state_i <- counties[i,'state_province']
    county_i <- counties[i,'county']
    
    centroid_i <-
      centroid_data |>
      dplyr::filter( country == country_i &
                state_province==state_i &
                county_parish==county_i )
    
    occurrences_i <-
      occurrences |>
      dplyr::filter( country == country_i &
                state_province==state_i &
                county==county_i )
    
    sf::st_crs(occurrences_i) <- st_crs(NULL)
    sf::st_crs(centroid_i) <- st_crs(NULL)
    
    if(nrow(centroid_i) > 0){
      
      
      
      
      # dists_m <- st_distance(occurrences_i, centroid_i)
      # 
      # dists |>
      #   apply(MARGIN = 1, FUN = which.min) -> dist_m_index
      # 
      # dists |>
      #   apply(MARGIN = 1, FUN = min) -> mins_m
      
      
      #set the crs to null so that the units will be decimal degrees (to match the max distances used in BIEN)
      
      sf::st_crs(occurrences_i) <- sf::st_crs(NULL)
      sf::st_crs(centroid_i) <- sf::st_crs(NULL)
      
      #calculate distances  
      
      dists_dd <- sf::st_distance(occurrences_i, centroid_i)
      
      dists_dd |>
        apply(MARGIN = 1, FUN = which.min) -> dist_dd_index
      
      dists_dd |>
        apply(MARGIN = 1, FUN = min) -> mins_dd
      
      closest_centroids_i <- centroid_i[dist_dd_index,]
      
      closest_centroids_i <-
        closest_centroids_i |>
        dplyr::mutate(cent_ID = occurrences_i$cent_ID,
               centroid_distance = mins_dd,
               relative_distance = centroid_distance/dist_max)|>
        sf::st_drop_geometry()
      
      if(nrow(closest_centroids_i) != nrow(occurrences_i)){
        stop("occurrence vs centroid mismatch in counties")}
      
      
      county_centroids |>
        dplyr::bind_rows(closest_centroids_i) -> county_centroids
      
      
    }else{
      
      #stop("Do something if there aren't any centroids (NAs prob)")
      
      
      occurrences_i |>
        dplyr::select(cent_ID, country, state_province, county) |>
        dplyr::rename(county_parish = county) |>
        sf::st_drop_geometry() |>
        dplyr::bind_rows(county_centroids) -> county_centroids
      
      
    }
    
    
  }#i county loop
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #return if no country is found 
  # if (nrow(country_centroids)==0) {
  #   
  # }
    
  
  
  # combining output
  country_centroids |>
    dplyr::bind_rows(state_centroids) |>
    dplyr::bind_rows(county_centroids) -> all_centroids
  
  # add in NA's where no political information was included
  if (nrow(all_centroids) == 0) stop ('no political information matched country centroids')
  
  occurrences |> 
    sf::st_drop_geometry() |>
    dplyr::filter(!cent_ID %in% (all_centroids |> dplyr::pull(cent_ID)))|>
    dplyr::select(cent_ID)|>
    dplyr::bind_rows(all_centroids) -> all_centroids
  
  
  # check lengths
  
  if(nrow(occurrences) != nrow(all_centroids)){stop("wrong number of centroids")}
  
  
  #reorder and rename fields

  #@BM REMOVED centroid_type of the list bc not in the initial centroid db
    all_centroids |>
      dplyr::relocate(c("cent_ID","id","gid_0","gid_1","gid_2","country",
                        "state_province","county_parish",#"centroid_type",
                        "relative_distance","centroid_distance","dist_max")) |>
      dplyr::arrange(cent_ID) |>
      tibble::remove_rownames() |>
      dplyr::rename(county = county_parish)->    all_centroids
  
  # columns we don't need and that might cause join issues
  
    all_centroids |>
      dplyr::select("cent_ID","relative_distance",#"centroid_type",
             "centroid_distance","dist_max") -> all_centroids
  
  # combine with occurrence data

    dplyr::full_join(occurrences,
              all_centroids,
              dplyr::join_by(cent_ID)) -> occurrences

    
  # evaluate relative distance vs threshold
    
    occurrences |>
      dplyr::mutate(is_centroid = dplyr::case_when(relative_distance <= relative_distance_threshold ~ 1,
                                     relative_distance > relative_distance_threshold ~ 0)) -> occurrences
  
  # select only the important columns
    
    occurrences |>
      dplyr::select("cent_ID","is_centroid","relative_distance",#"centroid_type",
             "centroid_distance","dist_max")|>
      dplyr::select(-cent_ID)|>
      sf::st_drop_geometry(.data) -> occurrences
  
  return(occurrences)
  
} #end fx
