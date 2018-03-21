\name{plotinitial}
\alias{plotinitial}
\title{
Initial Values plot.
}
\description{
Plot the fited curve and data, at a certain point of parameters, a given user parameter value (initial values mostly), or from \code{selfStart} slot.
}
\usage{
plotinitial(form, data, start = getInitial(form, data), length.out = 100,...)
}
\arguments{
  \item{form}{
\code{\link{nl.form}} object of the nonlinear regression model function.
}
  \item{data}{
List of data set.
}
  \item{start}{
NULL or list of starting value, include the parameter names.
}
  \item{length.out}{
Length of the incremented predictor values, to acheive a more smooth curve.
}
  \item{...}{extra argument to be passed to plot function, for further developement purposes.}
}
\details{
Plot initial can plot any fited curve at any parameter value. Technically it is implemented to test initial values computing by \code{selfStart} slot of \code{\link{nl.form}} object. But user can insert hi/her own starting value or even the final fited value.
}
\value{
Two dimensional Fitted curve and points.
%% ...
}
\references{
Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.

Seber, G., A. F. and Wild, C. J. (2003). Nonlinear Regression. New York: John Wiley & Sons, Inc.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
Although the final parameter estimates can be used to fit the curve, but it has less options than \code{plot} function.
}
\seealso{
 \code{\link{nl.form}}
 }
\examples{
## 
## Plot initial values from selfStart of Scaled Exponential Convex for carbon data
##
  crbdt<-list(xr=nlr::carbon$year,yr=nlr::carbon$co2)
  	plotinitial(form = nlrobj5[[8]],data = crbdt)
}
\keyword{ Initial Values }

