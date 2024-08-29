### Miscellanous functions  basic functions ###

# .subsetlist.nonNULL ====
#' @title Subset to non-NULL lists
#' @description Subset a list removing NULL entries
#' @param x List
#' @return list
#' @keywords internal
#' @author Josep M Serra-Diaz 
#' @noRd
.subsetlist.nonNULL <- function (x) {base::Filter(Negate(is.null), x)}

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
    stop("Issue: no output columns *_test exist. Maybe occTest filter results initially due to quality issues?")
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
