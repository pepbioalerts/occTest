% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Checks.R
\name{.addmainfields2}
\alias{.addmainfields2}
\title{Add main fields}
\usage{
.addmainfields2(dat, species, verbose = FALSE)
}
\arguments{
\item{dat}{A dataframe containing occurrence data for checking.}

\item{species}{character. Name of the species}

\item{verbose}{logical. Print messages? Default FALSE}
}
\value{
Original dataframe, dat.
}
\description{
Incorporate fields in the initial data frame
}
\details{
Inspired biogeo-package addmainfields but modified (hence number 2 after the function original name)
}
\seealso{
Other checks: 
\code{\link{.check.geospatial.data}()},
\code{\link{.checkdatastr2}()},
\code{\link{.checkfields}()},
\code{\link{.status.tracker.and.escaping}()}
}
\author{
Mark Robertson and Vernon Visser (original function), Josep M Serra Diaz (modifs)
}
\concept{checks}
