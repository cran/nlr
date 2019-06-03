\name{loss.hetroWM}
\alias{loss.hetroWM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Weighted Robut loss function.
}
\description{
  weighted loss function is used to estimate Weighted M-estimate. It is a robustified form of Likelihood function for heteroscedastic variance case. \code{loss.hetroWM} will be used in optimization function for estimating parameters of nonlinear function model and variance function model simultaneously. 
}
\usage{
loss.hetroWM(formula, data, start, varmodel, robfunc, ...)
}
\arguments{
  \item{formula}{
    \code{\link{nl.form}} object of nonlinear regression model.
}
  \item{data}{
list of data include responce, predictor or possibly predictor of variance model function (t), if not represented then the predict of nonlinear model function will be replaced in predictor variable of nonlinear variance model function that is \eqn{Var(\varepsilon)=\sigma^2 H(f(\theta),\tau)}}
  \item{start}{
    list of parameter values of nonlinear model function (\eqn{\theta} in \eqn{f(x,\theta)} and \eqn{\tau} in \eqn{H(x,\tau)}), initial values or increament during optimization procedure. 
}
  \item{varmodel}{
\code{\link{nl.form}} object of heteroscedastic variance model function. \code{varmodel} must return variance function not standard deviation. Include gradient and hessian.
}
  \item{robfunc}{
    \code{\link{nl.form}} of rho function. It must include tuning constants k0 and k1.
}
  \item{\dots}{
    extra argument might pass to nonlinear regression, heteroscedastic functions, or tuning constant and other parameters to robust loss rho function.
}
}
\details{
  Loss function in general form have robust rho function. 
  l(tau,theta)=sum[log(h(mu;tau,sg))] + sum[rho(ri/h(mu;tau,sg))]
}
\value{
  result <- list (value=value,angvec=angvec,angmat=angmat,
		refvar=refvar,sourcefnc= match.call(),
			rho=robvalue,fmod=fmod,varcomp=varcomp,correlation =nlrho,ri=ri)
  \item{value }{Robustified form of log-likelihood function, use in optimize function. Include gradient and hessian attributes.}
  \item{angvec }{angular vector for checking the convergence.}
  \item{angmat }{angular matrix for checking convergence in optimization procedure.}
  \item{refvar }{refvar, sample variance \eqn{ \sum(wi * z) \over \sum(wi)}}
  \item{sourcefnc}{Object of class \code{"callorNULL"} source function called for fitt.}
  \item{rho}{computed robust loss rho function. Include hessian and gradients.}
  \item{fmod}{computed function contains esponse and or its gradient and hessian predictor and or its gradient & hessian}
  \item{varcomp}{computed variance function model \eqn{H(t,\tau)}, include response or predictor and their gradient and hessian depends on the defined form of nonlinear function model.}
  \item{correlation}{correlation of fitt.}
  \item{ri}{residuals}
}
\references{
Lim, C., Sen, P. K., Peddada, S. D. (2010). Statistical inference in nonlinear regression under heteroscedasticity. Sankhya B 72:202-218.
}
\author{
Lim, C., Sen, P. K., Peddada, S. D. 
}
\note{
rho function can be square function, this will produce non robust MLE or LS.
This is implemented for internal use, might not be called directly by user.
}
\seealso{
  \code{\link{nl.form}}
}
\examples{
## The function is currently defined as
"loss.hetroWM"
}
\keyword{Wighted M-estimate}
