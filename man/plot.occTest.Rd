% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.occTest}
\alias{plot.occTest}
\title{Display the filtering process}
\usage{
\method{plot}{occTest}(x, occFilter_list = NULL, show_plot = FALSE, ...)
}
\arguments{
\item{x}{An \emph{occTest} object returned by  {\link[=occTest]{occTest}}, i.e. the unfiltered data.frame}

\item{occFilter_list}{Optional. An  \emph{occFilter} object; a list returned by  {\link[=occFilter]{occFilter}}, the result of the filtering of \code{x}}

\item{show_plot}{\emph{Logical} , should the plots be plotted ?}

\item{...}{not used}
}
\value{
\emph{list} of ggplots objects, of varying length, depending on whether the filtering was done by testBlock or testType
}
\description{
Display the filtering process
}
\details{
If \code{occFilter_list} is provided, display how the occurrences passed the different tests, otherwise only plot the coordinates filtering step
}
\examples{
#load output from occTest
occTest_output <- readRDS (system.file('ext/output_occTest.rds',package = 'occTest'))
#filter dataset output from occTest
filtered_occTest <- occFilter (occTest_output)
#plot the outputs
descriptive_plots <- plot (x=occTest_output,occFilter_list=filtered_occTest)
}
\seealso{
{\link[=occFilter]{occFilter}}  , {\link[=occTest]{occTest}}  , the {\link[=ggplot2]{ggplot2}} package
}
\author{
Jeremy Borderieux (jeremy.borderieux@agroparistech.fr)
}
\keyword{filter}
\keyword{plot}
