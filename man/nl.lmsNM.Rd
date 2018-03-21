\name{nl.lmsNM}
\alias{nl.lmsNM}
\title{
  Fitt a nonlinear regression model by least median of squares. The Optimization is done by golden section method.
}
\description{
  Least Median of square estimate is robust fitt by minimizing the median of squared residuals. This function use the "optim" function which minimize using Golden section method.
}
\usage{
nl.lmsNM(formula, data, start)
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
}
\details{
    Note that due to using classic optimization method, and due to non uniqueness of minimum of objective function this estimate is not efficient, but is high breakdown estimate. 
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
  A more efficient function is recomended is \code{\link{nl.lmsGA}}. When that function is not working properly this function can be used, but requires more precise initial values. These function uses by \code{\link{nlr}} for initial purposes. User can request least median square as initial value from \code{\link{nlr.control}} argument by \code{initials="lms"} method, nlr use \code{\link{nl.lmsGA}} not \code{nl.lmsNM}.
}

\seealso{
  \code{\link{nl.lmsGA}},  \code{\link{nlr.control}},  \code{\link{nlr}}
}
\examples{
  # chicken data fitt example
  data=list(xr=Weights$Date,yr=Weights$Weight)
  fit<- nl.lmsNM(nlrobj1[[14]],data=data,start=list(p1=1000,p2=42,p3=.11))
  fit
}
\keyword{ least median square estimate }
\keyword{ ordinary least square }
\keyword{ genetic algorithm }