% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{fullSettings}
\alias{fullSettings}
\title{Load full settings for occTest}
\usage{
fullSettings()
}
\value{
a \emph{list} with the parameter settings to run \link[occTest]{occTest}
}
\description{
Loads a list of lists with the different default parameters for analysi. 
It avoids using  some functions of the pkg under development.
}
\details{
it can be use internally or it can be used by a user to subsequently modify parameters
}
\examples{
\dontrun{
#load default settings
settings <- fullSettings()
}
}
\author{
Josep M Serra-Diaz (pep.serradiaz@agroparistech.fr)
}
\keyword{user}
