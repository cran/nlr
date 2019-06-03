\name{loss.SSQ}
\alias{loss.SSQ}
\title{ 
  Sum of squared loss function.
}
\description{
  This function used in \code{\link{nlsnm}} function to compute the least square estimate using derivative free Nelder-Mead algorithm. 
}
\usage{
loss.SSQ(formula, data, start, vm = NULL, rm = NULL, ...)
}

\arguments{
  \item{formula}{
        \code{\link{nl.form}} object of nonlinear regression model.
}
  \item{data}{
    list of data include responce and predictor.
}
  \item{start}{
     list of parameter values of nonlinear model function (\eqn{\theta} in \eqn{f(x,\theta)}), initial values or increament during optimization procedure.}
  \item{vm}{
    optional covariance matrix.
}
  \item{rm}{
    optional cholesky decomposition of covariance matrix.
}
  \item{\dots}{
        any other arguments might be used in formula, robfunc or tuning constants in rho function.
}
}
\details{
    \code{loss.SSQ} compute the sum of square of residuals, it is optimized to be used in \code{\link{nlsnm}} function, since optimization method Nelder-Mead is derivative free the result does not include derivatives.
}
\value{
  result <- list(value = value,correlation=correlation,fmod=fmod)

list values:
  \item{value }{sum of squared residuals.}
  \item{correlation }{correlation of model}
  \item{fmod}{computed function (transformed by R) contains esponse and or its gradient and hessian predictor and or its gradient & hessian, transformed also by R.}
}
\references{
  Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
 If required to compute square loss function include can use \code{nl.robfuncs[7]}, see \code{\link{nl.robfuncs}}.
 This is implemented for internal use, might not be called directly by user.
}
\seealso{
   \code{\link{nlsnm}}, \code{\link{nl.robfuncs}}
}
\examples{
## The function is currently defined as
"loss.SSQ"
}

\keyword{Nonlinear least square estimate}
\keyword{Nelder-Mead algorithm}
