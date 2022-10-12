#' @import ggplot2
#' @import sf
#' @import stringr
NULL

#' @title Display the filtering process
#' @descriptions Display the filtering process 
#' @details If \code{occFilter_list} is provided, display how the occurences passed the different tests, otherwise only plot the coordinates filtering step
#' @return list of ggplots objects, of varying length, depending on whether the filtering was done by testBlock or testType
#' @keywords plot filter
#' @author Jeremy Borderieux (jeremy.borderieux@@agroparistech.fr)
#' @param x An  occTest object returned by  {\link[=occTest]{occTest}}, i.e. the unfiltered data.frame 
#' @param occFilter_list Optional, an occFilter object; a list returned by  {\link[=occFilter]{occFilter}}, the result of the filtering of \code{x}
#' @param show_plot  Logical, should the plots be plotted ?
#' @param ... not used
#' @seealso  {\link[=occFilter]{occFilter}}  , {\link[=occTest]{occTest}}  , the {\link[=ggplot2]{ggplot2}} package
#' @examples 
#' #load output from occTest
#' occTest_output <- readRDS (system.file('ext/out.rds',package = 'occTest'))
#' #filter dataset output from occTest
#' filtered_occTest <- occFilter (occTest_output)
#' #plot the outputs
#' descriptive_plots <- plot (x=occTest_output,occFilter_list=filtered_occTest)
#' @export plot.occTest
#' @export

plot.occTest<-function(x,occFilter_list=NULL,show_plot=FALSE,...){
  ## check for S2 usage within sf, disabling it for the call if TRUE
  s2_used<-FALSE
  if(sf::sf_use_s2()){
    sf::sf_use_s2(FALSE)
    s2_used<-TRUE
  }
  
  ## get settings
  Settings<- get_occTest_settings(x)
  tableSettings<-Settings$tableSettings
  analysisSettings<-Settings$analysisSettings
  x_field<-tableSettings$x.field
  y_field<-tableSettings$y.field
  points_CRS<-analysisSettings$geoSettings$points.proj4string
  
  
  full_dataset<-x
  
  n_coords_missing<-sum(full_dataset$coordIssues_coordMissing_value)
  n_coords_filtered<-sum(full_dataset$Exclude)-n_coords_missing
  
  # we extract up to date country boundaries, but if it fails, we have a local copy of the shape_file (2020)
  countries_natural_earth<-try(st_as_sf(rnaturalearth::ne_countries(scale=50)),silent = TRUE)
  if("try-error" %in% class(countries_natural_earth)) {
    dest_url = 'https://github.com/pepbioalerts/vignetteXTRA-occTest/raw/main/ext/Pays_WGS84.rds'
    outFile = paste0(tempdir(),'/Pays_WGS84.rds')
    if (!file.exists(outFile)) utils::download.file(url=dest_url,destfile = outFile)
    countries_natural_earth<-readRDS(outFile)
  }
  
  
  if(raster::compareCRS(points_CRS,countries_natural_earth))countries_natural_earth<-st_transform(countries_natural_earth,crs=st_crs(points_CRS))
  
  ## cran check return notes for variable used with the ggplot syntaxe used within a function
  ## we can disable si note by declaring these vairable beforehand
  ## https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Exclude<-Reason<-occ_state<-NULL
  
  ## this dataframe contains ecery occurences, useful for the fist map displaying the first filtering
  ## we check for missing coordinated also
  full_dataset<-full_dataset[!full_dataset$coordIssues_coordMissing_value,]
  full_dataset$Reason<-ifelse(is.na(full_dataset$Reason),"Occurences kept for later tests",full_dataset$Reason)
  
  full_dataset_sf<-st_as_sf(full_dataset,coords=c(x_field,y_field),crs=st_crs(points_CRS))
  exten_of_plot_1<-st_bbox(full_dataset_sf)
  
  pch_filtered<-ifelse(n_coords_filtered>5000,17,4)# to fasten the display for bigger dataset, we replace crosses with triangles
  
  list_of_ggplot<-list()
  ## first plot: filtering by bad coordinates
  
  reasons<-unique(full_dataset[full_dataset$Reason!="Occurences kept for later tests",]$Reason)
  df_col<-data.frame(res=c("Occurences kept for later tests",reasons),
                     colors=c("grey45",c("brown","goldenrod","dodgerblue3","darkorchid3","forestgreen","lightblue3")[1:length(reasons)]),stringsAsFactors = FALSE)
  
  list_of_ggplot[[1]]<-ggplot(full_dataset_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
    geom_sf(aes(color=Reason,shape=Exclude==1),size=0.75,stroke=1.75,alpha=0.65)+
    theme(legend.position = "bottom")+
    coord_sf(xlim=c(exten_of_plot_1[1],exten_of_plot_1[3]),ylim=c(exten_of_plot_1[2],exten_of_plot_1[4]),expand = TRUE)+
    labs(color="Filter reason",title="Coordinates filtering",subtitle =paste0( "Occurences without coordinates = ",n_coords_missing ,"\nOccurences filtered =  ",n_coords_filtered))+
    scale_shape_manual(breaks = c(TRUE,FALSE),values=c(pch_filtered,16))+guides(shape = "none")+
    scale_color_manual(breaks =df_col$res ,values=df_col$colors)+ guides(colour = guide_legend(ncol = 2, byrow = TRUE,override.aes=list(pch=16)))
  
  if(is.null(occFilter_list)) warning("No occFilter object provided, plotting only the filtering by coordinates")
  if(!is.null(occFilter_list)){
    
    if (is.list(occFilter_list)) filtered_dataset<-occFilter_list$filteredDataset
    if (is.data.frame (occFilter_list) ) filtered_dataset<-occFilter_list
    if(nrow(filtered_dataset)==nrow(full_dataset)) message ("occFilter filtered 0 occurences")
    
    summary_stats<-occFilter_list$summaryStats
    
    filtered_dataset_scores<-cbind(full_dataset[full_dataset$Exclude==0,c(x_field ,y_field,"Exclude")],summary_stats)#"name"
    
    # we extract the score to run the filtering again
    dfScoreVals = filtered_dataset_scores %>% dplyr::select(dplyr::ends_with('_score')) 
    errorRule = occFilter_list$rule
    
    # @Pep avoiding CRAN-check notes on binding variables due to the use of tidyverse
    # filter_occ<-function(score,testName){
    #   current_treshold<-errorRule %>% dplyr::filter (test == stringr::str_remove(testName,"_score")) %>% dplyr::pull('errorThreshold')
    #   return( ifelse(is.na(score),FALSE, score>=current_treshold ))
    # }
    
    filter_occ<-function(score,testName){
      current_treshold<-errorRule 
      keepRows = which (current_treshold[,'test'] == stringr::str_remove(testName,"_score"))
      current_treshold = current_treshold [keepRows,]
      current_treshold = current_treshold %>% dplyr::pull('errorThreshold')
      return( ifelse(is.na(score),FALSE, score>current_treshold ))
      
    }
    
    toss_df<-mapply(filter_occ,dfScoreVals,colnames(dfScoreVals),SIMPLIFY = TRUE)
    colnames(toss_df)<-paste0(colnames(toss_df),"_bool")
    toss<-rowSums(toss_df)
    toss<-toss>=1
    
    
    filtered_dataset_scores<-cbind(filtered_dataset_scores,toss_df);rm(toss_df)
    filtered_dataset_scores$filtered<-toss
    
    ## This dataframe (transformed to an sf object) contain the filtered occurences, and which tests they passed or not
    filtered_dataset_scores_sf<-st_as_sf(filtered_dataset_scores,coords=c(x_field,y_field),crs=st_crs(points_CRS))
    # we will zomm in
    exten_of_plot_2<-st_bbox(filtered_dataset_scores_sf)
    
    
    ## we create this function, it will plot the result of every testBlock or testType passed
    plot_a_test_results<-function(score_name){
      
      ## avoid the note
      occ_state<-NULL
      
      
      ## test blocks get a special title
      title_char<-switch(score_name,geo_score="Geographical outliers tests",env_score="Environmental outliers tests",
                         lu_score="Land uses and human influence tests",time_score="Time precision test",sub("_score","",score_name))# Specific headers for the testBlocks
      
      bool_score_name<-paste0(score_name,"_bool")
      
      pch_filtered<-ifelse(sum(filtered_dataset_scores[,bool_score_name]>5000,na.rm=TRUE),17,4)
      
      
      if(all(is.na(filtered_dataset_scores[,bool_score_name]))) return(paste0("No '",sub("_score","",score_name),"' tests done"))
      filtered_dataset_scores_sf$occ_state<-ifelse(filtered_dataset_scores[,bool_score_name],"Removed",ifelse(filtered_dataset_scores$filtered,"Removed by another test","Kept occurences"))
      n_occ_filtered<-sum(filtered_dataset_scores[,bool_score_name])
      
      
      plot_result<-ggplot(filtered_dataset_scores_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
        geom_sf(aes(color=occ_state,shape=occ_state=="Removed"),size=0.75,stroke=1.75,alpha=0.5)+
        theme(legend.position = "bottom")+
        coord_sf(xlim=c(exten_of_plot_2[1],exten_of_plot_2[3]),ylim=c(exten_of_plot_2[2],exten_of_plot_2[4]),expand = TRUE)+
        labs(color="Filter results",title=title_char,subtitle =ifelse(n_occ_filtered==0,"No occurences filtered by these tests",paste0(n_occ_filtered," occurences removed by these tests")))+
        scale_shape_manual(breaks = c(TRUE,FALSE),labels=c("Rejected","Passed"),values=c(pch_filtered,16))+guides(shape = "none",colour=guide_legend(override.aes=list(pch=16)))+
        scale_color_manual(breaks =c("Removed","Removed by another test","Kept occurences"),values=c("brown","goldenrod","grey45") )
      
      return(plot_result)
    }
    
    test_executed<-colnames(dfScoreVals)
    
    ## we ran the function for each tests
    list_of_ggplot<-c(list_of_ggplot,lapply(test_executed,plot_a_test_results))
  }
  ## display to the user
  user_press<-"no_press"
  
  if(show_plot) for (current_plot in 1:length(list_of_ggplot)){
    print(list_of_ggplot[[current_plot]])
    if(current_plot!=length(list_of_ggplot))user_press<-readline(prompt="Press [enter] to see the next plot, type 'no' to stop plotting ")
    if(trimws(tolower(user_press))=="no")break
  }
  
  if(s2_used)sf::sf_use_s2(TRUE)
  
  return(list_of_ggplot)
  
}


#' @title Get occTest Settings
#' @description Get the settings used to create a  occTest or occFilter object
#' @return list of lists with all different parameters to use in  {\link[=occTest]{occTest}} 
#' @keywords occTest 
#' @author Jeremy Borderieux (jeremy.borderieux@@agroparistech.fr)
#' @param x An occTest or occFilter object returned by  {\link[=occTest]{occTest}} or  {\link[=occFilter]{occFilter}}
#' @seealso {\link[=occTest]{occTest}}; {\link[=occFilter]{occFilter}}
#' @examples 
#' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#'
#' #load output from occTest
#' occTest_output <- readRDS (system.file('ext/out.rds',package = 'occTest'))
#' get_occTest_settings(occTest_output)
#' @export 

get_occTest_settings<-function(x){
  if(!"occTest"%in%class(x) & !"occFilter"%in%class(x) ) warning("Not an occTest or an occFilter object, returning NULL instead of the list Settings")
  if(!"occTest"%in%class(x) & !"occFilter"%in%class(x) ) NULL  else  return(attr(x,"Settings"))
  
}
  
  

