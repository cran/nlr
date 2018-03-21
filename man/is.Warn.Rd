\name{is.Warn}
\alias{is.Warn}
\title{
  Check warning
}
\description{
  Return TRUE if warning occured. Error or warning of a procedure in \code{nlr} package will be saved in \code{Fault} slot of fited object with object type class \code{\link{Fault}}. \code{is.Warn} read the \code{Fault} slot or a \code{Fault} object to detect warning happened or no. The value inside \code{Fault$FN} will display the error code and \code{Fault$FL} represent error, both together can be used to identify warning.
}
\usage{
  is.Warn(obj)
}
\arguments{
  \item{obj}{
    all fited objects include a slot of \code{\link{Fault}} object to represent the error. \code{obj} can be fitt objects: \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{Fault}}.
  }
}
\details{
  \code{is.Warn} is imitating other \code{is.} primitive functions.It is implemented for debuging purpose in \code{nlr}. Warning might happen for example when maximum number of iteration exceeded, in this case the result might not be reliable and options should be changed. Another example might be when an infinit, null or rmissing value happened, in some cases not any result might be returned. For example if a derivative of a nonlinear function does not exist, the gradient values might have null values then the Nelder-Mead derivative free should be used.
}
\value{
  logical value \code{TRUE} or \code{T} if warning occured, \code{FALSE} or \code{F} if result is not warning. Technicaly if The slot \code{FN} of \code{Fault} object if is not zero means error or warning raised and \code{FL} is true if error occured. Other slots represent the error or warning informations.
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
  The \code{Fault} object can be error or warning, \code{\link{is.Faultwarn}} check for error or warning and \code{\link{is.Fault}} check the error.
}
\seealso{
  \code{\link{Fault}}, \code{\link{is.Faultwarn}} or \code{\link{is.Fault}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
is.Warn(Fault(FN=2))   # Not warning
is.Warn(Fault(FN=3))   # warning
}
\keyword{Debuging}
\keyword{Fault}
