# occFilter ====
#' @title Filter occurrence records from occTest outputs
#' @description Select occurrence records based on aggregated values of different tests
#' @return list of 2 data.frames
#' @keywords filter
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr), Jeremy Borderieux (jeremy.borderieux@@agroparistech.fr)
#' @param df data.frame. Output of  occTest
#' @param by character. Applying thresholds to either  blocks of test ('testBlock') or single test types ('testType')
#' @param errorAcceptance  character. Philosophy for filtering based on threshold. Option are majority, relaxed, strict. Default are 'relaxed'
#' @param errorThreshold double. Value from 0 to 1, specifying the threshold of wrong tests (potentally erroneous records) to filter. It overrides the parameters in thresholds. We recommend building that table based on the functio buildCustomThresholds.
#' @param custom data.frame or equivalent, custom rules created adding a "errorThreshold" (ranging from 0, strict, to 1, relaxed) column to to the result of readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
#' @details If errorAcceptance is used, a 'relaxed' philosophy corresponds to 0.7 (70% of tests of a block or type not passed), 'majority' corresponds to an errorAcceptance of 0.5, 'strict' corresponds to an errorAcceptance of 0.2.
#' @seealso showTests
#' @examples 
#' 
#' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#'
#' #load output from occTest
#' occTest_output <- readRDS (system.file('ext/out.rds',package = 'occTest'))
#' filtered_dataset <- occFilter (occTest_output)
#' #inspect results
#' names (filtered_dataset)
#' 
#' @export

occFilter <- function (df,
                       by='testBlock', #other option is testType
                       errorAcceptance = 'relaxed',
                       errorThreshold = NULL,
                       custom=NULL) {
  
  #load column metadata
  if (is.null(custom)){
    colMetaData = readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
    colMetaData = dplyr::as_tibble(colMetaData)
  }
  
  if (!is.null(custom)){
    colMetaData = dplyr::as_tibble(custom)
  }
  
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
  dfFiltered  = dplyr::as_tibble(df) 
  keep=which(dfFiltered[,'Exclude']==0)
  dfFiltered = dfFiltered[keep,]
             
        
  #get the scores for each category
  dfScores = lapply (uniqueTargetLevels, function (categ){
    baseCategSelec = myCategories %>% dplyr::filter (targetLevel==categ) %>% dplyr::pull(var = 'baseLevel')
    
    idCols = grep (pattern = paste0(baseCategSelec,collapse = '|'), names(dfFiltered),value = TRUE)
    idCols = grep (pattern ='_test', idCols,value = TRUE)
    
    dfSelec = dfFiltered [,idCols]
    
    dfScore= data.frame  (score= .gimme.score (dfSelec),
                          rateTestPerformed = apply (dfSelec,MARGIN = 1,function  (x) sum (!is.na(x))/length(x)),
                          NtestsPerformed =   apply (dfSelec,MARGIN = 1, function  (x) (sum (!is.na(x)))))
    names (dfScore) = paste0(categ,'_',names(dfScore))
    dfScore
  })
  dfScores = dplyr::bind_cols(dfScores)
  dfScoreVals = dfScores %>% dplyr::select(dplyr::ends_with('_score')) 
  
  #rules of selection 
  if (is.null(errorThreshold) & is.null(custom)){
    if (! errorAcceptance %in% c('strict','majority','relaxed')) stop (paste0('errorAcceptance',errorAcceptance,' type not known'))
    if (errorAcceptance == 'strict')    {colMetaData$errorThreshold =0.2}
    if (errorAcceptance == 'majority')  {colMetaData$errorThreshold =0.5}
    if (errorAcceptance == 'relaxed')   {colMetaData$errorThreshold =0.7}
   
  }
  if (is.null(errorAcceptance) & is.null(custom)){
    colMetaData$errorThreshold =errorThreshold
   
  }
  
  nDfScore = unlist(strsplit (names (dfScoreVals),split = '_score'))
   
  errorThresholdDf = colMetaData %>% 
      dplyr::filter (mode != 'filter') %>%
      dplyr::select(by,errorThreshold) %>%
      unique
   
  names (errorThresholdDf) <- c('test','errorThreshold')
  
  # errorThresholdDf = errorThresholdDf %>%
  # dplyr::filter(test %in% nDfScore)
  errorThresholdDf = as.data.frame (errorThresholdDf)
  keep_errorThresholdDf = which (errorThresholdDf[,'test'] %in% nDfScore)
  errorThresholdDf = errorThresholdDf [keep_errorThresholdDf,]
  
  if (!is.null(custom)){

   errorThresholdDf = errorThresholdDf[match(nDfScore, errorThresholdDf$test),]

   errorThreshold = errorThresholdDf %>% dplyr::pull (errorThreshold)

    if (length(errorThreshold) != ncol (dfScoreVals)) {stop (paste('Different threshold for the same',by))}
  
  }
    errorRule=errorThresholdDf
    
    # filter_occ<-function(score,testName){
    #   current_treshold<-errorRule %>% dplyr::filter (test == stringr::str_remove(testName,"_score")) %>% dplyr::pull(errorThreshold)
    #   return( ifelse(is.na(score),FALSE, score>current_treshold ))
    # 
    # }
    filter_occ<-function(score,testName){
      
      # current_treshold<-errorRule %>% dplyr::filter (test == stringr::str_remove(testName,"_score")) %>% dplyr::pull(errorThreshold)
      # return( ifelse(is.na(score),FALSE, score>current_treshold ))
      # 
      current_treshold<-errorRule 
      keepRows = which (current_treshold[,'test'] == stringr::str_remove(testName,"_score"))
      current_treshold = current_treshold [keepRows,]
      current_treshold = current_treshold %>% dplyr::pull('errorThreshold')
      return( ifelse(is.na(score),FALSE, score>current_treshold ))

    }


      toss_df<-mapply(filter_occ,
                      dfScoreVals,
                      colnames(dfScoreVals),SIMPLIFY = TRUE)
    
    # @JB I added na.rm=TRUE in case some scores are not done for a given record
    if (any (class(toss_df) %in% 'logical')) {toss= sum(toss_df,na.rm=TRUE)} else
    {toss<-rowSums(toss_df,na.rm=TRUE)}
      
    toss<-toss>=1
    
  #output
  if (all(toss==FALSE)) {warning("No occurence filtered through quality tests")
    out = list (filteredDataset = dfFiltered,
                summaryStats = dfScores,
                rule = errorRule)
    attr(out,"class")<-c("occFilter",class(out))
    attr(out,"Settings")<-get_occTest_settings(df)
    return(out)
    }
  if (all(toss==TRUE)) {
    warning("All occurence filtered")
    out = list (filteredDataset = NULL,
                summaryStats = dfScores,
                rule = errorRule)
    attr(out,"class")<-c("occFilter",class(out))
    attr(out,"Settings")<-get_occTest_settings(df)
    return(out)
    }
  out = list (filteredDataset = dfFiltered[which(toss)*(-1),],
              summaryStats = dfScores,
              rule = errorRule)
  attr(out,"class")<-c("occFilter",class(out))
  attr(out,"Settings")<-get_occTest_settings(df)
  return(out)
  
}



