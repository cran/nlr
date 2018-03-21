\name{is.Fault}
\alias{is.Fault}
\title{
  Check error
}
\description{
  Error of a procedure in \code{nlr} package will be saved in \code{Fault} slot of fited object with object type class \code{\link{Fault}}. \code{is.Fault} read the \code{Fault} slot or a \code{Fault} object to detect error happened or no. The value inside \code{Fault$FL} will show the result.
}
\usage{
is.Fault(obj)
}
\arguments{
  \item{obj}{
    all fited objects include a slot of \code{\link{Fault}} object to represent the error. \code{obj} can be fitt objects: \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{Fault}}.
  }
}
\details{
  \code{is.Fault} is imitating other \code{is.} primitive functions.It is implemented for debuging purpose in \code{nlr}.
}
\value{
  logical value \code{TRUE} or \code{T} if error happened, \code{FALSE} or \code{F} if result is not error.
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
  The \code{Fault} object can be warning, but not error, in this case you can use \code{\link{is.Warn}} or \code{\link{is.Faultwarn}} to detect warnings.
}
\seealso{
  \code{\link{Fault}}, \code{\link{is.Warn}} or \code{\link{is.Faultwarn}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
is.Fault(1)
is.Fault(Fault(FN=1))
}
\keyword{Fault}
\keyword{Debuging}
