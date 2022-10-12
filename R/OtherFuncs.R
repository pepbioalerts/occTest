### Miscellanous functions  basic functions ###



# .subsetlist.nonNULL ====
#' @title Subset to non-NULL lists
#' @description Subset a list removing NULL entries
#' @details
#' @param x List
#' @return list
#' @keywords internal
#' @author Josep M Serra-Diaz 
#' @noRd
.subsetlist.nonNULL <- function (x) {base::Filter(Negate(is.null), x)}



# .subsetlist.isNA ====
#' @title Subset to non-NA lists
#' @description Subset a list removing NA entries
#' @details
#' @param x List
#' @return list
#' @keywords internal
#' @author Josep M Serra Diaz
#' @noRd
.subsetlist.isNA <- function (x) {Filter(Negate(is.na), x)}





# .join.spname ====
#' @title Convert to Collated species names
#' @description It collates species names from format "Species genus" to "Species_genus"
#' @param x character
#' @return character
#' @keywords internal
#' @noRd
.join.spname <- function (x){
  a <- strsplit(x,split=' ')
  a <- unlist (a)
  a <- paste (a,collapse = '_')
  return(a)
}


# .multiple.strsplit  ====
#' @title multiple splitting of strings
#' @description uses several splitting characters to split a character string.
#' @param x character.
#' @param multiple.splits character. Vector with the characters used to recursively split x
#' @return character
#' @keywords internal
#' @author Pep Serra-Diaz
#' @noRd
# Function written by Andrew Bevan, found on R-sig-Geo, and modified by Pascal Title
.multiple.strsplit <- function (x,multiple.splits) {

  stopifnot( length(multiple.splits)>1 )

  a <- unlist( strsplit (x, split = multiple.splits[1] ) )

  for (i in 2:length(multiple.splits)){
    a <- unlist( strsplit (a, split = multiple.splits[i] ) )
  }

  return(a)
}

# .paste3  ====
#' @title Paste obviating NAs
#' @description Removes NAs when pasting character vectors
#' @param ... character vector
#' @return character
#' @keywords internal
#' @author Pep Serra-Diaz
#' @noRd

.paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

# .gimme.score  ====
#' @title get the scoring outputs from an analysis
#' @description compute the scores for each type of test implemented in occTest function
#' @details This is implemented within the occTest function
#' @param x data.frame resulting from a series of analysis in occTest
#' @return character
#' @keywords internal 
#' @author Pep Serra-Diaz
#' @family internal
#' @noRd
.gimme.score <- function (x){
  if (!is.data.frame(x)) {stop ('input needs to be a dataframe')}

  
  #identify the needed columns in the data frame
  columns.for.scoring = grep ('_test',names(x))
  
  if (length(columns.for.scoring)>1) {
    score <- rowMeans (x [,columns.for.scoring], na.rm = TRUE)
    
  }
  if (length(columns.for.scoring)==1){
    score<- x[,columns.for.scoring]
  }
  
  if (length(columns.for.scoring)==0){
    stop("Issue: seems that no output columns *_test exist")
  }
  
  score
  
  
}


# .get_os  ====
#' @title Identify type of Operative system
#' @description The function is used to report the working operative system as a 3-letter character for 'win', and 'mac', 'unix', and 'unknown OS'
#' @details This is implemented whan parallelizing functions
#' @return character
#' @keywords internal 
#' @author Pep Serra-Diaz
#' @note Internal function copied from utils wickham to know which is the os used by the computer
#' @seealso 
#' @aliases
#' @family internal

.get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}


# .find_and_transform_dates  ====
#' @title Identify and harmonize data formats
#' @description haromnize different date styles 
#' @return character
#' @keywords internal 
#' @seealso dataPreparation::find_and_transform_dates
#' @import dataPreparation 

.find_and_transform_dates <- function (data_set, cols = "auto", formats = NULL, n_test = 30, 
          ambiguities = "IGNORE", verbose = TRUE) {
  function_name <- "find_and_transform_dates"
  data_set <- check_and_return_datatable(data_set)
  is.verbose(verbose)
  is_ambiguities(ambiguities, function_name)
  cols <- real_cols(data_set = data_set, cols = cols, function_name = function_name)
  start_time <- proc.time()
  formats_found <- dataPreparation::identify_dates(data_set, cols = cols, formats = formats, 
                                  n_test = n_test, ambiguities = ambiguities, verbose = verbose)
  if (verbose) {
    printl(function_name, ": It took me ", round((proc.time() - 
                                                    start_time)[[3]], 2), "s to identify formats")
  }
  if (length(formats_found) < 1) {
    if (verbose) {
      printl(function_name, ": There are no dates to transform.\n                   (If i missed something please provide the date format in inputs or\n                   consider using set_col_as_date to transform it).")
    }
    return(data_set)
  }
  start_time <- proc.time()
  data_set <- dataPreparation::set_col_as_date(data_set, format = formats_found, 
                              verbose = FALSE)
  if (verbose) {
    printl(function_name, ": It took me ", round((proc.time() - 
                                                    start_time)[[3]], 2), "s to transform ", length(formats_found), 
           " columns to a Date format.")
  }
  return(data_set)
}


# splitSpname  ====
#' @title Split joined species name
#' @description  By default, it converts "Species_genus" format to "Species genus" format
#' @param x species name join as a string separated by "_" 
#' @return a species name where genus and species are separated by a space 
#' @family internal
splitSpname = function (x){
  a = strsplit(x,'_')[[1]]
  paste(a,collapse = ' ')
  
}

# rm.all =====
#' @title Remove all objects
#' @description  removes all objects in the environment
#' @return a clear environment
#' @family internal
rm.all <- function () {rm (list =  setdiff(ls(),'rm.all') ) ; gc(verbose = FALSE) }

# lazylogic =====
#' @title a lazy logic omitting NAs
#' @description  evaluates a conditional expression where if any element is NA the output is FALSE instead of NA
#' @param e expression 
#' @return logic
#' @family internal

lazylogic <- function (e){
  
  o <- eval(parse(text=e))
  ifelse(is.na(o),FALSE,o)
  
}

# .hijack [deprecated] ====
# #' @title Hijack functions
# #' @description  Hijacking functions to rename them
# #' @details In the occTest package this is used to get to a same function name for different OS implementing differnt parallelization systems
# #' @return character
# #' @keywords internal 
# #' @author tylerrinker on Rbloggers
# #' @note got it from https://www.r-bloggers.com/2014/08/hijacking-r-functions-changing-default-arguments/
# #' @family internal
# .hijack <- function (FUN, ...) {
#   .FUN <- FUN
#   args <- list(...)
#   invisible(lapply(seq_along(args), function(i) {
#     formals(.FUN)[[names(args)[i]]] <<- args[[i]]
#   }))
#   .FUN
#  }
