### SHOW FUNCTIONS
# show_tableSettings ====
#' @title Print naming conventions in occTest
#' @description prints a table with the the conventions used for column names
#' @details The function prints a guide to column naming conventions used by occTest in their default parameters. 
#'          These defaults can be changed via set_tableNames, but the user may also decide to format their input table according to these naming conventions. 
#'          It does not require input parameters
#' @returns prints a data.frame with the tableSettings parameters
#' @keywords user
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples 
#' ts <- show_tableSettings ()
#' @export
show_tableSettings <- function (){
  tabNames=readRDS(system.file('ext/tableColumns.rds',package='occTest'))
  DT::datatable(tabNames)
}



# show_writeoutSettings ====
#' @title Print naming conventions and definitions for writoutSettings list
#' @description prints a table with the the conventions used for writoutSettings list
#' @details The function prints a guide to the list naming conventions used by occTest in their default parameters. 
#'          These defaults can be changed via set_writeout. 
#'          It does not require input parameters
#' @return prints a data.frame with the writeout parameters
#' @keywords user
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples 
#' show_writeoutSettings ()
#' @export
show_writeoutSettings <- function (){
  tabNames=readRDS(system.file('ext/writeoutSettings_metadata.rds',package='occTest'))
  DT::datatable(tabNames)
}



# show_analysisSettings ====
#' @title Print naming conventions and definitions for analysisSettings list
#' @description prints a table with the the conventions used for analysisSettings list
#' @details The function prints a guide to the list naming conventions used by occTest in their default parameters. 
#'          These defaults can be changed manually by the use. To follow the same structure we recommend loading defaultSettings() and then modify $analysisSettings based on that structure 
#'          It does not require input parameters
#' @return prints a data.frame with the analysis settings
#' @keywords user
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples 
#' show_analysisSettings ()
#' @export
show_analysisSettings <- function (){
  tabNames=readRDS(system.file('ext/analysisSettings_metadata.rds',package='occTest'))
  DT::datatable(tabNames)
}


# show_tests ====
#' @title Show implemented tests and types of tests 
#' @description prints a table with the column names
#' @details The function prints a guide to column naming conventions used by occTest in their default parameters. These defaults can be changed via set_tableNames, but the user may also decide to format their input table according to these naming conventions. 
#' @return prints a dataframe with the avaialble tests
#' @keywords user
#' @author Josep M Serra-Diaz (pep.serradiaz@@agroparistech.fr)
#' @examples
#' show_tests()
#' @export
show_tests<- function (){
  tabNames=readRDS(system.file('ext/fieldMetadata.rds',package='occTest'))
  DT::datatable(tabNames)
}




