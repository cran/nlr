\name{is.Faultwarn}
\alias{is.Faultwarn}
\title{
  Check error or warning
}
\description{
  Return TRUE if error oe warning occured. Error or warning of a procedure in \code{nlr} package will be saved in \code{Fault} slot of fited object with object type class \code{\link{Fault}}. \code{is.Faultwarn} read the \code{Fault} slot or a \code{Fault} object to detect error happened or no. The value inside \code{Fault$FN} will display the error code. If is not zero means an error or warning raised.
}
\usage{
is.Faultwarn(obj)
}
\arguments{
  \item{obj}{
    all fited objects include a slot of \code{\link{Fault}} object to represent the error. \code{obj} can be fitt objects: \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{Fault}}.
  }
}
\details{
  \code{is.Faultwarn} is imitating other \code{is.} primitive functions.It is implemented for debuging purpose in \code{nlr}. 
}
\value{
  logical value \code{TRUE} or \code{T} if error or warning occured, \code{FALSE} or \code{F} if result is not error nor warning. Technicaly if The slot \code{FN} of \code{Fault} object if is not zero means error or warning raised, in that case other slots represent the error or warning informations. 
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  The \code{Fault} object can be error or warning, \code{\link{is.Warn}} check for warning only and \code{\link{is.Fault}} check the error.
}
\seealso{
  \code{\link{Fault}}, \code{\link{is.Warn}} or \code{\link{is.Fault}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
is.Faultwarn(Fault(FN=1))
}
\keyword{Debuging}
\keyword{Fault}
