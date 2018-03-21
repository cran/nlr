\name{nl.lmsGA}
\alias{nl.lmsGA}
\title{
  Fitt a nonlinear regression model by least median of squares. The Optimization is done by genetic algorithm.
}
\description{
  Least Median of square estimate is robust fitt by minimizing the median of squared residuals. This function use the "ga" function,from "GA" package, which minimize using genetic algorithm method.
}
\usage{
nl.lmsGA(formula, data, start,min=NULL,max=NULL,type="real-valued")
}

\arguments{
  \item{formula}{
      \code{\link{nl.form}} object of the nonlinear function model.
}
  \item{data}{
      list of data include responce and predictor.
}
  \item{start}{
      list of parameter values of nonlinear model function (\eqn{\theta}. in \eqn{f(x,\theta)}).
}
  \item{min}{vector of minimum values of parameters, which passes to "ga" function.}
  \item{max}{vector of maximum values of parameters, which passes to "ga" function.}
  \item{type}{the type of genetic algorithm to be run, which passes to "ga" function.}
}
\details{
    Note that due to using genetic algorithm method, this function is mor efficient than \code{\link{nl.lmsNM}}, but estimator due to non uniqueness of minimum of objective function is not efficient, but is high breakdown estimate. 
}
\value{
  list of parameter estimates.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  When that function is not working properly, \code{\link{nl.lmsNM}} function can be used, but requires more precise initial values. These function uses by \code{\link{nlr}} for initial purposes. User can request least median square as initial value from \code{\link{nlr.control}} argument by \code{initials="lms"} method, nlr use \code{\link{nl.lmsGA}} not \code{nl.lmsNM}.
}

\seealso{
  \code{\link{nl.lmsNM}},  \code{\link{nlr.control}},  \code{\link{nlr}}
}
\examples{
  # chicken data fitt example
  data=list(xr=Weights$Date,yr=Weights$Weight)
  fit<- nl.lmsGA(nlrobj1[[14]],data=data,start=list(p1=100,p2=42,p3=.11))
}
\keyword{ least median square estimate }
\keyword{ ordinary least square }
\keyword{ genetic algorithm }