\name{dfr.robloss}
\alias{dfr.robloss}
\title{
  Robut loss function.}
\description{
  Resturn robust loss function for minimization purpose to find the M-estimate. It is used in \code{\link{dfrmest.NLM}} function for derivative free purpose. Gradient and hessian are computed numerically. 
}
\usage{
dfr.robloss(formula, data, start, robfunc, control = nlr.control(), rmat = NULL, ...)
}
\arguments{
  \item{formula}{
    \code{\link{nl.form}} object of nonlinear regression model.
}
  \item{data}{
    list of data include responce and predictor.
}
  \item{start}{
    list of parameter values of nonlinear model function (\eqn{\theta} in \eqn{f(x,\theta)}), initial values or increament during optimization procedure. It must include scale sigma (standard deviation), if not included Fault(9) will be returned.
}
  \item{robfunc}{
    \code{\link{nl.form}} of rho function. It must include tuning constants k0 and k1.
}
  \item{control}{
    list of \code{\link{nlr.control}} for controling convergence criterions.
}
  \item{rmat}{R-Matrix for transforming, it might be cholesky decomposition of covariance matrix.}
  \item{\dots}{
    any other arguments might be used in formula, robfunc or tuning constants in rho function.
}
}
\details{
  Compute Loss function, sum of robust rho function to compute the M-estimate.
  \deqn{\ell(\theta)=\sum \rho\left(\frac{r_i}{\sigma}\right)}
  Standard deviation \eqn{\sigma} must be included in \code{start} argument list with the name \code{sigma}.
  
  \code{gradient} and \code{hessian} attributes compute numerically.
}
\value{  result <- list(htheta=htheta,rho=robvalue,ri=rsd,fmod=fmod,Fault=Fault2)
  list of output:
  \item{htheta:}{sum of rho function, include attribute \code{"gradient"} and \code{"hessian"}}
  \item{rho:}{computed rho function and attributes of \code{"gradient"} and \code{"hessian"}}
  \item{ri:}{residuals}
  \item{fmod:}{computed function contains esponse and or its gradient and hessian predictor and or its gradient & hessian}
  \item{Fault:}{\code{\link{Fault}} object of error, if no error Fault number = 0 will return back.}
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
All functions should have \code{gradient} and \code{hessian} in attributes. For derivative free purpose the \code{\link{dfr.robloss}} can be used.
It is designed for internal use, might not call directly by user.
}
\seealso{
\code{\link{nl.form}}, \code{\link{nlr.control}}, \code{\link{nlmest.NLM}}, \code{\link{dfr.robloss}}, \code{\link{robloss}}
}
\examples{
## The function is currently defined as
"dfr.robloss"
}
\keyword{M-estimate}
\keyword{Robust Loss function}
\keyword{Optimization}
\keyword{Location estimate}
\keyword{Nelder-Mead optimization}
