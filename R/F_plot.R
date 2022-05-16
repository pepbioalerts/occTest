#### display ####

#' @title Display the filtering process
#' @descriptions Select occurrence records based on aggregated values of different tests
#' @details
#' @param 
#' @return list of 2 data.frames
#' @keywords plot filter
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @param df data.frame outputed from occurrenceClassify
#' @param by character. Applying thresholds to either  blocks of test ('testBlock') or single test types (testTypes )
#' @param threshold  character. Phylosophy for filtering based on threshold. Option are majority, relaxed, stringent 
#' @param customThresholds data.frame. Specify for each dimension of flags the specific threshold to use. It overides the parameters in thresholds. We recomend building that table based on the functio buildCustomThresholds.
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
#' @export

plot_occFilter<-function(occfilter_list,occTest_df,show_plot=T){
  
  
  filtered_dataset<-occfilter_list$filteredDataset
  
  if(nrow(filtered_dataset)==nrow(occTest_df))stop("occFilter filtered 0 occurences")
  
  summary_stats<-occfilter_list$summaryStats
  full_dataset<-occTest_df
  treshold<-as.numeric(occfilter_list$rule[2])
  
  filtered_dataset_scores<-cbind(full_dataset[full_dataset$Exclude==0,c("decimalLongitude" ,"decimalLatitude","Exclude")],summary_stats)#"name"
  
  n_coords_missing<-sum(full_dataset$coordIssues_coordMissing_value)
  n_coords_filtered<-sum(full_dataset$Exclude)-n_coords_missing
  
  dfScoreVals = filtered_dataset_scores %>% dplyr::select(ends_with('_score')) 
  filtered_dataset_scores$filtered<- apply (dfScoreVals,1,function (x) {
    a = any (x >= treshold)
    if (is.na(a)) a <- F
    a
  })
  
  
  
  tmp_score_passed<-apply (dfScoreVals,2,function (x)x<treshold)
  colnames(tmp_score_passed)<-paste0(colnames(tmp_score_passed),"_bool")
  
  filtered_dataset_scores<-cbind(filtered_dataset_scores,tmp_score_passed);rm(tmp_score_passed)
  filtered_dataset_scores_sf<-st_as_sf(filtered_dataset_scores,coords=c("decimalLongitude","decimalLatitude"),crs=st_crs(4326))
  
  countries_natural_earth<-try(st_as_sf(rnaturalearth::ne_countries(scale=50)),silent = T)
  if("try-error" %in% class(countries_natural_earth))countries_natural_earth<-read_sf(system.file('ext/Country_shp_world',package='occTest'),layer="Pays_WGS84")
  
  
  exten_of_plot_2<-st_bbox(filtered_dataset_scores_sf)
  
  
  list_of_ggplot<-list()
  
  
  
  full_dataset_sf<-st_as_sf(full_dataset[!full_dataset$coordIssues_coordMissing_value,],coords=c("decimalLongitude","decimalLatitude"),crs=st_crs(4326))
  full_dataset_sf$Reason<-ifelse(is.na(full_dataset_sf$Reason),"Occurences kept for later tests",full_dataset_sf$Reason)
  exten_of_plot_1<-st_bbox(full_dataset_sf)
  
  list_of_ggplot[[1]]<-ggplot(full_dataset_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
    geom_sf(aes(color=Reason,shape=Exclude==1),size=0.75,stroke=1.75,alpha=0.65)+
    theme(legend.position = "bottom")+
    coord_sf(xlim=c(exten_of_plot_1[1],exten_of_plot_1[3]),ylim=c(exten_of_plot_1[2],exten_of_plot_1[4]),expand = T)+
    labs(color="Filter reason",title="Coordinates filtering",subtitle =paste0( "Occurences without coordinates = ",n_coords_missing ,"\nOccurences filtered =  ",n_coords_filtered))+
    scale_shape_manual(breaks = c(TRUE,FALSE),values=c(4,16))+guides(shape = FALSE)+
    scale_color_manual(values=c("brown","goldenrod","dodgerblue3","grey45"))+ guides(colour = guide_legend(ncol = 2, byrow = T))
  
  
  
  
  
  
  plot_a_test_results<-function(score_name){
    title_char<-switch(score_name,geo_score="Geographical outliers tests",env_score="Environental outliers tests",
                       lu_score="Land uses and human influence tests",time_score="Time precision test",sub("_score","",score_name))# Specific headers for the testBlocks
    
    bool_score_name<-paste0(score_name,"_bool")
    
    if(all(is.na(filtered_dataset_scores[,bool_score_name]))) return(paste0("No '",sub("_score","",score_name),"' tests done"))
    filtered_dataset_scores_sf$occ_state<-ifelse(!filtered_dataset_scores[,bool_score_name],"Removed",ifelse(filtered_dataset_scores$filtered,"Removed by another test","Kept occurences"))
    n_occ_filtered<-sum(!filtered_dataset_scores[,bool_score_name])
    
    
    plot_result<-ggplot(filtered_dataset_scores_sf)+theme_bw()+geom_sf(data=countries_natural_earth,color="grey60",fill="grey98")+
      geom_sf(aes(color=occ_state,shape=occ_state=="Removed"),size=0.75,stroke=1.75,alpha=0.65)+
      theme(legend.position = "bottom")+
      coord_sf(xlim=c(exten_of_plot_2[1],exten_of_plot_2[3]),ylim=c(exten_of_plot_2[2],exten_of_plot_2[4]),expand = T)+
      labs(color="Filter results",title=title_char,subtitle =ifelse(n_occ_filtered==0,"No occurences filtered by these tests",paste0(n_occ_filtered," occurences removed by these tests")))+
      scale_shape_manual(breaks = c(TRUE,FALSE),labels=c("Rejected","Passed"),values=c(4,16))+guides(shape = FALSE)+
      scale_color_manual(breaks =c("Removed","Removed by another test","Kept occurences"),values=c("brown","goldenrod","grey45") )
    
    
    return(plot_result)
  }
  
  
  test_executed<-colnames(dfScoreVals)
  
  list_of_ggplot<-c(list_of_ggplot,lapply(test_executed,plot_a_test_results))
  
  
  
  if(show_plot) for (current_plot in 1:length(list_of_ggplot)){
    print(list_of_ggplot[[current_plot]])
    if(current_plot!=length(list_of_ggplot))user_press<-readline(prompt="Press [enter] to see the next plot, type 'no' to stop plotting ")
    if(trimws(tolower(user_press))=="no")break
  }
  
  
  return(list_of_ggplot)
  
  
}

