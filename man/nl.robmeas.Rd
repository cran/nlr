\name{nl.robmeas}
\Rdversion{1.1}
\docType{class}
\alias{nl.robmeas}
\alias{nl.robmeas-class}
\alias{$,nl.robmeas-method}
\alias{plot,nl.robmeas,ANY-method}
\alias{plot,nl.robmeas-method}
\title{Class \code{"nl.robmeas"}}
\description{
  Outlier detection measure object.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nl.robmeas", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{measure}:}{Object of class \code{"numeric"} vector of statistics measure for each data points.}
    
    \item{\code{cutofpoint}:}{Object of class \code{"numeric"} cut of point for the measure.}
    \item{\code{name}:}{Object of class \code{"character"} name of the measure.}
  }
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "nl.robmeas")}: access the slots.}
    \item{plot}{\code{signature(x = "nl.robmeas", y = "ANY")}: plot sequence of measure and cut of point line. Usually if the value is more than cut of point the data point might be outlier or influence observation.
      In contrast  \code{plot.atyps} function is to plot \code{"nl.robmeas"} object for multiple purpose. If more than such objects are stored in a list the plot will iterate on all cells and draw all. It can be called as:
      
      \code{plot.atyps(x,...)} #where x is lis that include more than one \code{"nl.robmeas"} object.
    }
	 }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  \code{nl.robmeas} usually generate from atypical method to calculate several measures for identifying outlier or influential observations.
}
\seealso{
  \code{\link{nlout}}, \code{\link{nl.fitt}}
}
\examples{
showClass("nl.robmeas")
}
\keyword{classes}
\keyword{outliers}
\keyword{influence}
\keyword{influential observatiosn}


