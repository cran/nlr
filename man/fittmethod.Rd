\name{fittmethod}
\Rdversion{1.1}
\docType{class}
\alias{fittmethod}
\alias{fittmethod-class}
\alias{$,fittmethod-method}
\title{Class \code{"fittmethod"}}
\description{
  Fitt method object, it store the method used for estimation, runing function.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fittmethod", ...)}.
}
\section{Slots}{
  \describe{
  \item{\code{methodID}:}{ Object of class \code{"numeric"} code for the method.}
  
   \item{\code{method}:}{Object of class \code{"character"} name of the method used in estimation and other computation procedures.}
   
   \item{\code{detail}:}{Object of class \code{"character"} detail text description of the method.}
   
  \item{\code{methodBR}:}{Object of class \code{"numeric"} (method Branche) branch for the method used in iteration.}
  
  \item{\code{detailBR}:}{Object of class \code{"character"} detail of the branch.}
  
  \item{\code{subroutine}:}{Object of class \code{"character"} the function, subroutine, that the result constructed from.}
  
    \item{\code{lossfunction}:}{Object of class \code{"character"} objective loss function used by subroutine to optimize.}
    
    \item{\code{subroutineBR}:}{Object of class \code{"character"} sub subroutine, called for optimize.}
  }
}
\section{Extends}{
Class \code{"\link{fittmethodorNULL}"}, directly.
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "fittmethod")}: access slots.}
	 }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  This is for internal use, might not be called directly by user. All the fited objects, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rgn}}, \code{\link{nl.fitt.rob}}, have a method slot of object type \code{fittmethod}. It can be used by \code{\link{recalc}} method to recalculate the fitt.
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rgn}}, \code{\link{nl.fitt.rob}}
}
\examples{
fittmethod(methodID=1)
}
\keyword{classes}
