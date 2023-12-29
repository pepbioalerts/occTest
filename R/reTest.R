### addTest =====
#'Add tests to an occTest output
#'
#'Reorganize an output from occTest columns to incorporate user defined tests.
#' @param occTest_result \emph{data.frame} or \emph{tibble}. Output from an occTest project
#' @param my_test \emph{data.frame} or tibble. New tests following the syntax of occTest (see Details)
#' @returns \emph{data.frame} with the new user-defined tests incorporated (see Details)
#' @details
#' my_test needs to follow the syntax of occTest testType_testMethod_test, testType_testMethod_value, testType_testMethod_comment 
#' the parameter my_test needs to have at least a logical column with the convention [testType]_[testMethod]_[test]
#'  where testTypes are the different tests types implemented in occTest
#' Note that the result has not have the testType_scores update. For that you need to use reTest.
#' @examples \donttest{
#' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#' occTest_output <- readRDS (system.file('ext/output_occTest.rds',package = 'occTest'))
#' myNewTest <- data.frame (humanDetection_myInventedTest_test=sample (c(TRUE,FALSE),
#' replace = TRUE,
#' size = nrow(occTest_output)))
#' out_test_newAdded <- addTest(occTest_result = occTest_output, my_test = myNewTest)
#' }
#' @export
addTest <- function (occTest_result, my_test){
  #identify columns in each dataset
  id_col_pattern <- function (ds,n=3){
    names_ds <- strsplit (names (ds),'_')
    names_ds <- sapply (names_ds, function (x) {
      if (length(x)==n) {x[1]} else {NA}
    })
    names_ds
  }
  names_occTest <- id_col_pattern (ds = occTest_result)
  names_new_tests <- id_col_pattern (ds = my_test)
  #check that all testTypes are in occTest 
  if (!all(names_new_tests %in% names_occTest)) stop ('currently addTest can only work on testTypes included in occTest')
  
  #identify start and end of occTest testTypes identified in new_test
  colsID_occTest <- lapply (unique(names_new_tests), function (x){which (names_occTest %in% x)})
  names (colsID_occTest) <-  unique(names_new_tests)
  
  #reorder columns in the sequence of occTest result
  colsID_occTest <- colsID_occTest[order (sapply (colsID_occTest,min))]
  
  #add colums
  for (i in 1:length(colsID_occTest)){
    if (i==1) out <- occTest_result
    new_testMethods <- grep (x = names(my_test) ,pattern = names (colsID_occTest)[i],value = T)
    #new_testMethods <- grep (x = new_testMethods , pattern = '_test$',value = T)
    out <- tibble::add_column(out,my_test[,new_testMethods,drop=F],.after=max (colsID_occTest[[i]]))
  }
  
  return (out)
}



### reTest =====
#'Updates occTest with new user-defined tests
#'
#'Reorganize an output from occTest columns to incorporate user defined tests.
#' @param occTest_result \emph{data.frame} or tibble. Output from an occTest project
#' @param my_new_test \emph{data.frame} or \emph{tibble}. New tests following the syntax of occTest (see Details)
#' @returns \emph{data.frame with} the new user-defined tests incorporated in the table and the testType_scores columns updated.
#' @details
#' my_test needs to follow the syntax of occTest testType_testMethod_test, 
#' testType_testMethod_value, testType_testMethod_comment 
#' the parameter my_test needs to have at least a logical column with 
#' the convention [testType]_[testMethod]_[test]
#' where testTypes are the different tests types implemented in occTest
#' @examples \donttest{
#' ### THIS IS A CUT DOWN  EXAMPLE 
#' ### visit vignetteXtra-occTest for more info
#' out <- readRDS (system.file('ext/output_occTest.rds',package = 'occTest'))
#' myNewTest <- data.frame (humanDetection_myInventedTest_test=sample (c(TRUE,FALSE),
#'                                                                     replace = TRUE,
#'                                                                     size = nrow(out)))
#' out_test_new <- reTest(occTest_result = out, my_new_test = myNewTest)
#' out_test_new
#' }
#' @export
reTest <- function (occTest_result, my_new_test){
  myNew_occTest <- addTest(occTest_result = occTest_result ,my_test = my_new_test)
  names_score <- grep (names (myNew_occTest),pattern='_score$',value = T)
  names_score <- unlist (strsplit (names_score,split = '_score'))
  for (n in names_score){
    idCols <- grep (names (myNew_occTest),pattern = paste0('^',n,'.*_test$'),value=T)
    myNewScore <- .gimme.score(myNew_occTest[,idCols,drop=F])
    myNewScore [which (is.nan(myNewScore))] <- NA
    myNew_occTest [, paste0(n,'_score')] <- myNewScore
  }
  return (myNew_occTest)
}

