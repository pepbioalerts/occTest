####### FILTERING OCCURRENCE 



#' @title Build tables of user specified filters 
#' @descriptions This function is used to build a table of customThresholds specified by the user to be used in occFilter
#' @details
#' @param 
#' @return data.frame
#' @keywords filter
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @param df data.frame outputed from occurrenceClassify
#' @param level numeric Level of application of the filtering process.
#' @note
#' @seealso
#' @references
#' @aliases
#' @family
#' @examples \dontrun{
#' example<-"goes here"
#' }
# #' @export

buildCustomThresholds <- function (df,level){
  
  #not implemented
  

  
  
}


#' @title Filter occurrence records
#' @descriptions Select occurrence records based on aggregated values of different tests
#' @details
#' @param 
#' @return list of 2 data.frames
#' @keywords filter
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @param df data.frame outputed from occurrenceClassify
#' @param level numeric. Level of application of the filtering process. 1 is the coarsest 3 is the minimum
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

occFilter <- function (df,
                       by='testBlock', #other option is testType
                       errorAcceptance = 'relaxed',
                       errorThreshold = NULL) {

   #load pipes (need to change)
   usethis::use_pipe()
   library (magrittr)

  #load column metadata
   #colMetaData = read.csv  (system.file('ext/fieldMetadata.csv',package='occTest'))
   colMetaData = readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
   colMetaData = dplyr::as_tibble(colMetaData)
   
   targetLevel = colMetaData %>% 
     dplyr::filter (mode != 'filter') %>%
     dplyr::select(by)
   
   basecateg = colMetaData %>% 
     dplyr::filter (mode != 'filter') %>%
     dplyr::select('testType')
   
   myCategories = dplyr::bind_cols(targetLevel,basecateg) %>% unique
   names (myCategories) = c('targetLevel','baseLevel')
   
   uniqueTargetLevels = unique (myCategories$targetLevel)
   
   #select target rows (not the ones with coord issues) 
   dfFiltered  = dplyr::as_tibble(df) %>%
     dplyr::filter (Exclude ==0 ) 
   
   #get the scores for each category
   dfScores = lapply (uniqueTargetLevels, function (categ){
     baseCategSelec = myCategories %>% dplyr::filter (targetLevel==categ) %>% dplyr::pull(var = baseLevel)
     
     dfSelec = dfFiltered %>% 
       dplyr::select(starts_with(baseCategSelec)) %>% 
       dplyr::select(ends_with('_test'))
       
       
      dfScore= data.frame  (score= .gimme.score (dfSelec),
               rateTestPerformed = apply (dfSelec,MARGIN = 1,function  (x) sum (!is.na(x))/length(x)),
               NtestsPerformed =  apply (dfSelec,MARGIN = 1,function  (x) (sum (!is.na(x)))))
       names (dfScore) = paste0(categ,'_',names(dfScore))
       dfScore
   })
  dfScores = dplyr::bind_cols(dfScores)
   
  dfScoreVals = dfScores %>% dplyr::select(ends_with('_score')) 
   
  #rules of selection 
  if (is.null(errorThreshold)){
    if (errorAcceptance == 'strict') {errorThreshold =0.2}
    if (errorAcceptance == 'relaxed')   {errorThreshold =0.7}
    if (errorAcceptance == 'majority')  {errorThreshold =0.5}
  }
  toss = apply (dfScoreVals,1,function (x) {any (x >= errorThreshold)})
  #output
  list (fitleredDataset = dfFiltered[!toss,],
        summaryStats = dfScores)
}
