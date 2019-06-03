\name{robloss}
\alias{robloss}
\title{
  Robut loss function.
}
\description{
Resturn robust loss function for minimization purpose to find the M-estimate.
}
\usage{
robloss(formula, data, start, robfunc, control = nlr.control(robscale = T), ...)
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
  \item{\dots}{
    any other arguments might be used in formula, robfunc or tuning constants in rho function.
}
}
\details{
  Compute Loss function, sum of robust rho function to compute the M-estimate.
  \deqn{\ell(\theta)=\sum \rho\left(\frac{r_i}{\sigma}\right)}
  Standard deviation \eqn{\sigma} must be included in \code{start} argument list with the name \code{sigma}.
}
\value{
  list of output:
  \item{htheta}{sum of rho function, include attribute \code{"gradient"} and \code{"hessian"}}
  \item{rho}{computed rho function and attributes of \code{"gradient"} and \code{"hessian"}}
  \item{ri}{residuals}
  \item{hessh.p1}{hessian of loss function part1}
  \item{hessh.p2}{hessian of loss function part2, in clasic this part removed but in robust statistics values are significant and can not be omited, See Riazoshams et al 1014}
  \item{dtilda}{D(thilda) part of hessian}
  \item{fmod}{computed function contains esponse and or its gradient and hessian predictor and or its gradient & hessian}
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
  All functions should have \code{gradient} and \code{hessian} in attributes. For derivative free purpose the \code{\link{dfr.robloss}} can be used.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}

\seealso{
\code{\link{nl.form}}, \code{\link{nlr.control}}, \code{\link{nlmest.NLM}}, \code{\link{dfr.robloss}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"robloss"
}
\keyword{M-estimate}
\keyword{Robust Loss function}
\keyword{Optimization}
\keyword{Location estimate}
