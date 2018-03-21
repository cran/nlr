\name{loss.robchis}
\alias{loss.robchis}
\title{
  Hetroscedastic chi-square robust loss function.
}
\description{Compute the robustified chi-square pseoudo log likelihood,based on varaince model function.
\deqn{\ell (\boldsymbol{\theta},\tau )=\sum\{ w_i \log (H({x_i};\tau) + \rho \left[\sqrt {{z_i}/ H({x_i};{\sigma ^2},\lambda )} \right] \}}
}
\usage{
  loss.robchis(formula, data, start, theta, varmodel, robfunc, ...)
}

\arguments{
  \item{formula}{
    \code{\link{nl.form}} object of nonlinear function \eqn{f(x;\theta)}.
}
  \item{data}{
    list of data include responce, predictor or possibly predictor of variance model function (t), if not represented then the predict of nonlinear model function will be replaced in predictor variable of nonlinear variance model function that is \eqn{Var(\varepsilon)=\sigma^2 H(f(\theta),\tau)}
}
  \item{start}{
    list of parameter values of variance model function (\eqn{\tau} in \eqn{H(t,\tau)}), initial value or increament during optimization procedure.
}
  \item{theta}{
    list of model function parameter (\eqn{\theta} in \eqn{f(x,\theta)}).
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
    For estimating variance model parameter \eqn{\tau}, robustified form of chi-square pseudo likelihood is used as robust estimate. Based on calculating the sample variances.
}
\value{  
  list od los function values:
  \item{value}{value of minus loglikelihood of chi-square, include attribute "gradient"" and "hessian". These values use in optimization functions.}
  \item{angvec}{angular vector for checking the convergence.}
  \item{angmat}{angular matrix for checking convergence in optimization procedure.}
  \item{refvar}{refvar, sample variance \eqn{ \sum(wi * z) \over \sum(wi)}}
  \item{fmod}{computed function model \eqn{f(x,\theta)}, include response, predictor and their gradient and hessian depends on the defined form of nonlinear function model.}
  \item{varcomp}{computed variance function model \eqn{H(t,\tau)}, include response or predictor and their gradient and hessian depends on the defined form of nonlinear function model.}
  \item{vcmdata}{list of data used in variance model function, that is varmodel$independent and varmodel$dependent typically is zi.}
  \item{sourcefnc}{source function from which this function is called. May be used in feature computing such as outlier detection measures.}
  \item{rho}{computed robust loss rho function. Include hessian and gradients.}
  \item{zi}{computed sample variance, which follows the chi-square distribution.}
}
\references{
Riazoshams H 2010 Outlier detection and robust estimation methods for nonlinear regression having autocorrelated and heteroscedastic
errors.
}
\author{
 Hossein Riazoshams, 08/01/2010.
 Email: \email{riazihosein@gmail.com} 
 URL \url{hriaz.amarefars.com/nlr}
}
\note{
  This is used for robust estimates, for classic see \code{\link{loss.chis}}
  This is implemented for internal use, might not be called directly by user.  
}

\seealso{
  \code{\link{nl.form}}, \code{\link{loss.chis}}
}
\examples{
## The function is currently defined as
"loss.robchis"
}
\keyword{ heteroscedastic variance }
\keyword{ nonlinear regression function model}
\keyword{heteroscedastic function model}
\keyword{chi-square pseudo log likelihood}
\keyword{Robust estimators}

