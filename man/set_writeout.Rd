% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{set_writeout}
\alias{set_writeout}
\title{Set output objects from ocTest}
\usage{
set_writeout(
  output.dir = NULL,
  writeAllOutput = NULL,
  write.simple.output = NULL,
  write.full.output = NULL,
  output.base.filename = NULL
)
}
\arguments{
\item{output.dir}{character. Output directory}

\item{writeAllOutput}{logical. Should all outputs be written to a file? Overides write.simple.output and write.full.output}

\item{write.simple.output}{logical. Should only summary outputs be written to a file?}

\item{write.full.output}{logical. Should all detailed outputs be written to a file?}

\item{output.base.filename}{character. Name of the country code field.}
}
\value{
a \emph{list} with the parameter settings to run \link[occTest]{occTest}
}
\description{
helper function to set the options for the outputs (writeoutSettings) in occTest function. 
Defaults values may be found 
Alternatively, the user can specify their own field names for the table
}
\examples{
defaultOutput <- set_writeout()
#only modifying the names for the coordinates
myTable_withMyOutput <- set_writeout(writeAllOutput=TRUE) 
}
\author{
Josep M Serra-Diaz (pep.serradiaz@agroparistech.fr)
}
\keyword{user}
