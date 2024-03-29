% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Analysis.R
\name{duplicatesexcludeAnalysis}
\alias{duplicatesexcludeAnalysis}
\title{Duplicated records}
\usage{
duplicatesexcludeAnalysis(
  df,
  xf,
  yf,
  resolution.in.minutes,
  raster.grid = NULL,
  verbose = FALSE
)
}
\arguments{
\item{df}{\emph{data.frame} of species occurrences}

\item{xf}{the column name in df containing the x coordinators}

\item{yf}{the column name in df containing the y coordinators}

\item{resolution.in.minutes}{the resolution of environmental data used, specified in minutes}

\item{raster.grid}{An optional \emph{SpatRast} object}

\item{verbose}{\emph{logical}. Print messages? Default FALSE}
}
\value{
\emph{list} of three components: Dups.Exact = exact duplicate records, Dups.Grid= duplicates within environmental gridcell, continue = dataframe with good records

a \emph{list} of 3 \emph{data.frames} containing exact coordinte duplicates, grid duplicates (coordinates falling into the same grid cell) and non duplicated coordinates.
}
\description{
checks for duplicated coordinates in the occurrence data frame
}
\details{
it differentiates the exact duplicates and the duplicates for a occurrences falling in the same pixel
}
\examples{
ds <- data.frame (x=c(runif (n = 100),1000),y=c(runif (n = 100),1000),Reason=NA)
duplicatesexcludeAnalysis(ds,xf='x',yf='y',resolution.in.minutes=60)
}
\seealso{
Other analysis: 
\code{\link{.nearestcell3}()},
\code{\link{biogeo_outliers}()},
\code{\link{centroidDetection}()},
\code{\link{countryStatusRangeAnalysis}()},
\code{\link{geoEnvAccuracy}()},
\code{\link{humanDetection}()},
\code{\link{timeAccuracy}()}
}
\author{
Josep M Serra-Diaz (pep.serradiaz@agroparistech.fr)
}
\concept{analysis}
\keyword{internal}
