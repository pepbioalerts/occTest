# occFilter ====
#' @title Filter occurrence records from occTest outputs
#' @descriptions Select occurrence records based on aggregated values of different tests
#' @details
#' @param 
#' @return list of 2 data.frames
#' @keywords filter
#' @author JM Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @param df data.frame. Output of  occTest
#' @param by character. Applying thresholds to either  blocks of test ('testBlock') or single test types ('testType')
#' @param errorAcceptance  character. Philosophy for filtering based on threshold. Option are majority, relaxed, stringent. Default are 'relaxed'
#' @param errorThreshold double. Value from 0 to 1, specifying the threshold of wrong tests (potentally erroneous records) to filter. It overrides the parameters in thresholds. We recommend building that table based on the functio buildCustomThresholds.
#' @details If errorAcceptance is used, a 'relaxed' philosophy corresponds to 0.7 (70% of tests of a block or type not passed), 'majority' corresponds to an errorAcceptance of 0.5, 'stringent' corresponds to an errorAcceptance of 0.2.
#' @note
#' @seealso showTests
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

  #load pipes (need to change, you do not call libraries in functions)
  # usethis::use_pipe()
  # library (magrittr)

  #load column metadata
   colMetaData = readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
   colMetaData = dplyr::as_tibble(colMetaData)
   
   targetLevel = colMetaData %>% 
     dplyr::filter (mode != 'filter') %>%
     dplyr::select(by)
   
   basecateg = colMetaData %>% 
     dplyr::filter (mode != 'filter') %>%
     dplyr::select('method')
   
   myCategories = dplyr::bind_cols(targetLevel,basecateg) %>% unique
   names (myCategories) = c('targetLevel','baseLevel')
   
   uniqueTargetLevels = unique (myCategories$targetLevel)
   
   #select target rows (not the ones with coord issues) 
   dfFiltered  = dplyr::as_tibble(df) %>%
     dplyr::filter (Exclude ==0 ) 
   
   #get the scores for each category
   dfScores = lapply (uniqueTargetLevels, function (categ){
     baseCategSelec = myCategories %>% dplyr::filter (targetLevel==categ) %>% dplyr::pull(var = baseLevel)
     
     idCols = grep (pattern = paste0(baseCategSelec,collapse = '|'), names(dfFiltered),value = T)
     idCols = grep (pattern ='_test', idCols,value = T)
     
     dfSelec = dfFiltered [,idCols]
       
      dfScore= data.frame  (score= .gimme.score (dfSelec),
               rateTestPerformed = apply (dfSelec,MARGIN = 1,function  (x) sum (!is.na(x))/length(x)),
               NtestsPerformed =   apply (dfSelec,MARGIN = 1, function  (x) (sum (!is.na(x)))))
       names (dfScore) = paste0(categ,'_',names(dfScore))
       dfScore
   })
  dfScores = dplyr::bind_cols(dfScores)
  dfScoreVals = dfScores %>% dplyr::select(ends_with('_score')) 
   
  #rules of selection 
  if (is.null(errorThreshold)){
    if (! errorAcceptance %in% c('strict','majority','relaxed')) stop (paste0('errorAcceptance',errorAcceptance,' type not known'))
    if (errorAcceptance == 'strict')    {errorThreshold =0.2}
    if (errorAcceptance == 'majority')  {errorThreshold =0.5}
    if (errorAcceptance == 'relaxed')   {errorThreshold =0.7}
  }
  #select those rows with error rate above threshold, and ignore if output is NA (that means effectively it is a F)
  toss = apply (dfScoreVals,1,function (x) {
     a = any (x >= errorThreshold)
     if (is.na(a)) a <- F
     a
  })
  
  #output
  if (all(toss==F)) {
     out  = list (filteredDataset = dfFiltered,
           summaryStats = dfScores,
           rule = c(errorAcceptance,errorThreshold))
     attr(out,"class")<-c("occFilter",class(out))
     attr(out,"Settings")<-get_occTest_settings(df)
     
     return (out)
  } 
  
 out = list (filteredDataset = dfFiltered[which(toss)*(-1),],
        summaryStats = dfScores,
        rule = c(errorAcceptance,errorThreshold))
 attr(out,"class")<-c("occFilter",class(out))
 attr(out,"Settings")<-get_occTest_settings(df)
 return(out)
 
}

