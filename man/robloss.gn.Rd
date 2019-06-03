\name{robloss.gn}
\alias{robloss.gn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generalized Robut loss function.
}
\description{
Resturn Generalized robust loss function for minimization purpose to find the Generalized M-estimate. Generalized M-estimate required correlation or covariance matrix of data, then the model transform and estimated.
}
\usage{
robloss.gn(formula, data, start, robfunc, rmat, control = nlr.control(robscale = T), ...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{rmat}{
    R matrix, is cholesky decomposition of covariance matrix, the model transform by multiplying by R matrix.
}
  \item{control}{
    list of \code{\link{nlr.control}} for controling convergence criterions.
}
  \item{\dots}{
    any other arguments might be used in formula, robfunc or tuning constants in rho function.
}
}
\details{
 Compute Loss function, sum of robust rho function to compute the M-estimate.
  \deqn{\ell(\theta)=\sum \rho\left(\frac{R \times r_i}{\sigma}\right)}
  Standard deviation \eqn{\sigma} must be included in \code{start} argument list with the name \code{sigma}.
  
  The R matrix is \code{rmat} argument.
}
\value{
  list of output:
  \item{htheta}{sum of rho function, include attribute \code{"gradient"} and \code{"hessian"}}
  \item{rho}{computed rho function and attributes of \code{"gradient"} and \code{"hessian"}}
  \item{ri}{residuals, transformed by R.}
  \item{hessh.p1}{hessian of loss function part1}
  \item{hessh.p2}{hessian of loss function part2, in clasic this part removed but in robust statistics values are significant and can not be omited, See Riazoshams et al 1014}
  \item{dtilda}{D(thilda) part of hessian}
  \item{fmod}{computed function (transformed by R) contains response and or its gradient and hessian predictor, transformed also by R.}
  \item{Fault}{\code{\link{Fault}} object of error, if no error Fault number = 0 will return back.}
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{
  This function use in optimization functions, specially from \code{\link{nlmest.NLM}}, for where the covariance matrix or R matrix given, may not be called explicitly by user.
  
  Generalized M-estimate might reperesent the autocorrelated or heteroscedastic variance case.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
\code{\link{nl.form}}, \code{\link{nlr.control}}, \code{\link{nlmest.NLM}}
}
\examples{
## The function is currently defined as
"robloss.gn"
}
\keyword{M-estimate}
\keyword{Robust Loss function}
\keyword{Optimization}
\keyword{Location estimate}
\keyword{Heteroscedastic variance}
\keyword{Autocorrelated error}

