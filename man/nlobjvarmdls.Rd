\name{nlrobjvarmdls1}
\alias{nlrobjvarmdls1}
\alias{nlrobjvarmdls2}
\alias{nlrobjvarmdls3}
\docType{data}
\title{
  Variance model objects.
}
\description{
  List of \code{\link{nl.form}} objects of pre defined variance models. They use in heteroscedastic variance cases.
  \itemize{
    \item{nlrobjvarmdls1}{ has 7 nonlinear variance model function, with \eqn{\sigma^2} as variance.}
    \item{nlrobjvarmdls1}{ has 7 nonlinear variance model function exactly same as \code{nlrobjvarmdls1} but the \eqn{\sigma} as variance. In fact it can be standard deviation, square roots. But \code{nlr} work with variance functions.}
    \item{nlrobjvarmdls1}{ has 3 nonlinear variance model function, in general case they dont include constant variance \eqn{\sigma}. Variance is general parameteric form \eqn{Var(error)=H(x,\tau)}}
  }
}
\usage{
  data(nlrobjvarmdls1) 
  
  #nlrobjvarmdls1[[1]] access first element.
  }
\format{
The format is: chr "nlrobjvarmdls1"
 
 Elements are list of model, each cell is \code{\link{nl.form}} object.
}
\details{
  by index can access any of the nodels, for example nlrobj1[[14]] with name nlrobjvarmdls1[[1]]$name="power" is power model used to fitt the chicken \code{\link{Weights}} data.
}
\source{
Riazoshams 2015
}
\references{
  Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al.
}
\seealso{
  \code{\link{nl.form}}, \code{\link{Weights}}
}

\examples{
data(nlrobjvarmdls1)
nlrobjvarmdls1
}
\keyword{datasets}
\keyword{hetroscedastic variance}
