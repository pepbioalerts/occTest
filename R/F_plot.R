#' @import ggplot2
#' @import sf
NULL

#' @title Display the filtering process
#' @descriptions Display the filtering process 
#' @details If \code{occfilter_list} is provided, display how the occurences passed the different tests, otherwise only plot the coordinates filtering step
#' @return list of ggplots objects, of varying length, depending on whether the filtering was done by testBlock or testType
#' @keywords plot filter
#' @author J Borderieux (jeremy.borderieux@@agroparistech.fr)
#' @param occTest_df A data.frame returned by  {\link[=occTest]{occTest}}, i.e. the unfiltered data.frame 
#' @param occfilter_list Optional, a list returned by  {\link[=occFilter]{occTest}}, the result of the filtering by \code{occTest_df}
#' @param tableSettings  list, Elements corresponding to different settings of the input occurrence table, used to create the \code{occTest_df} object, if null, try to use the default settings
#' @param show_plot  Logical, should the plots be plotted ?
#' @note
#' @seealso  {\link[=occFilter]{occFilter}}  , {\link[=occTest]{occTest}}  , the {\link[=ggplot2]{ggplot2}} package
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export

plot_occFilter<-function(occTest_df,occfilter_list=NULL,tableSettings=NULL,show_plot=T){
  if(is.null(tableSettings))tableSettings<- defaultSettings()$tableSettings
  x_field<-tableSettings$x.field
  y_field<-tableSettings$y.field
  
  full_dataset<-occTest_df
  
  n_coords_missing<-sum(full_dataset$coordIssues_coordMissing_value)
  n_coords_filtered<-sum(full_dataset$Exclude)-n_coords_missing
  
  # we extract up to date contrie boundaries, but if it fails, we have a local copy of the shape_file (2020)
  countries_natural_earth<-try(st_as_sf(rnaturalearth::ne_countries(scale=50)),silent = T)
  if("try-error" %in% class(countries_natural_earth))countries_natural_earth<-read_sf(system.file('ext/Country_shp_world',package='occTest'),layer="Pays_WGS84")
  
  
  ## this dataframe contains ecery occurences, useful for the fist map displaying the first filtering
  ## we check for missing coordinated also
  full_dataset<-full_dataset[!full_dataset$coordIssues_coordMissing_value,]
  full_dataset$Reason<-ifelse(is.na(full_dataset$Reason),"Occurences kept for later tests",full_dataset$Reason)
  
  full_dataset_sf<-st_as_sf(full_dataset,coords=c(x_field,y_field),crs=st_crs(4326))## automate the crs selection
  exten_of_plot_1<-st_bbox(full_dataset_sf)
  
  
  list_of_ggplot<-list()
  ## first plot: filtering by bad coordinates
  
  reasons<-unique(full_dataset[full_dataset$Reason!="Occurences kept for later tests",]$Reason)
  df_col<-data.frame(res=c("Occurences kept for later tests",reasons),
                     colors=c("grey45",c("brown","goldenrod","dodgerblue3","darkorchid3","forestgreen","lightblue3")[1:length(reasons)]),stringsAsFactors = F)
  
  list_of_ggplot[[1]]<-ggplot(full_dataset_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
    geom_sf(aes(color=Reason,shape=Exclude==1),size=0.75,stroke=1.75,alpha=0.65)+
    theme(legend.position = "bottom")+
    coord_sf(xlim=c(exten_of_plot_1[1],exten_of_plot_1[3]),ylim=c(exten_of_plot_1[2],exten_of_plot_1[4]),expand = T)+
    labs(color="Filter reason",title="Coordinates filtering",subtitle =paste0( "Occurences without coordinates = ",n_coords_missing ,"\nOccurences filtered =  ",n_coords_filtered))+
    scale_shape_manual(breaks = c(TRUE,FALSE),values=c(4,16))+guides(shape = FALSE)+
    scale_color_manual(breaks =df_col$res ,values=df_col$colors)+ guides(colour = guide_legend(ncol = 2, byrow = T))
  
  
  if(!is.null(occfilter_list)){
  
  filtered_dataset<-occfilter_list$filteredDataset
  
  if(nrow(filtered_dataset)==nrow(occTest_df))stop("occFilter filtered 0 occurences")
  
  summary_stats<-occfilter_list$summaryStats
  treshold<-as.numeric(occfilter_list$rule[2])
  
  filtered_dataset_scores<-cbind(full_dataset[full_dataset$Exclude==0,c(x_field ,y_field,"Exclude")],summary_stats)#"name"

  # we extract the score to run the filtering again
  dfScoreVals = filtered_dataset_scores %>% dplyr::select(ends_with('_score')) 
  filtered_dataset_scores$filtered<- apply (dfScoreVals,1,function (x) {
    a = any (x >= treshold)
    if (is.na(a)) a <- F
    a
  })
  
  
  # we run the filtering again, and transform each tests passed as a logical
  tmp_score_passed<-apply (dfScoreVals,2,function (x)x<treshold)
  colnames(tmp_score_passed)<-paste0(colnames(tmp_score_passed),"_bool")
  
  filtered_dataset_scores<-cbind(filtered_dataset_scores,tmp_score_passed);rm(tmp_score_passed)
 

  ## This dataframe (transformed to an sf object) contain the filtered occurences, and which tests they passed or not
  filtered_dataset_scores_sf<-st_as_sf(filtered_dataset_scores,coords=c(x_field,y_field),crs=st_crs(4326))
  # we will zomm in
  exten_of_plot_2<-st_bbox(filtered_dataset_scores_sf)
  
  
  ## we create this function, it will plot the result of every testBlock or testType passed
  plot_a_test_results<-function(score_name){
    ## test blocks get a special title
    title_char<-switch(score_name,geo_score="Geographical outliers tests",env_score="Environental outliers tests",
                       lu_score="Land uses and human influence tests",time_score="Time precision test",sub("_score","",score_name))# Specific headers for the testBlocks
    
    bool_score_name<-paste0(score_name,"_bool")
    
    if(all(is.na(filtered_dataset_scores[,bool_score_name]))) return(paste0("No '",sub("_score","",score_name),"' tests done"))
    filtered_dataset_scores_sf$occ_state<-ifelse(!filtered_dataset_scores[,bool_score_name],"Removed",ifelse(filtered_dataset_scores$filtered,"Removed by another test","Kept occurences"))
    n_occ_filtered<-sum(!filtered_dataset_scores[,bool_score_name])
    
    
    plot_result<-ggplot(filtered_dataset_scores_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
      geom_sf(aes(color=occ_state,shape=occ_state=="Removed"),size=0.75,stroke=1.75,alpha=0.5)+
      theme(legend.position = "bottom")+
      coord_sf(xlim=c(exten_of_plot_2[1],exten_of_plot_2[3]),ylim=c(exten_of_plot_2[2],exten_of_plot_2[4]),expand = T)+
      labs(color="Filter results",title=title_char,subtitle =ifelse(n_occ_filtered==0,"No occurences filtered by these tests",paste0(n_occ_filtered," occurences removed by these tests")))+
      scale_shape_manual(breaks = c(TRUE,FALSE),labels=c("Rejected","Passed"),values=c(4,16))+guides(shape = FALSE)+
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
  
  return(list_of_ggplot)
  
}

