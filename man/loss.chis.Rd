\name{loss.chis}
\alias{loss.chis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Hetroscedastic chi-square loss function.
}
\description{
  Compute the minus of chi-square pseoudo log likelihood,based on varaince model function.
  \eqn{\ell ({\sigma ^2},\lambda ) = \sum \{ w_i \log (\tilde H({x_i};{\sigma ^2},\lambda )) + {z_i}/\tilde H({x_i};{\sigma ^2},\lambda ) \}}
}
\usage{
  loss.chis(formula, data, start, theta, varmodel, ...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{\dots}{
    extra argument might pass to nonlinear regression or heteroscedastic functions.
}
}
\details{
    For estimating variance model parameter \eqn{\tau}, chi-square pseudo chi square is used as classic estimate. Based on calculating the sample variances.
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
  \item{zi}{computed sample variance, which follows the chi-square distribution.}
}
\references{
  Bunke, O., Droge, B., Polzehl, J. Splus tools for model selection in nonlinear regression (1998) Computational Statistics, 13 (2), pp. 257-281.
}
\author{
  Bunke, O., Droge, B., Polzehl
}
\note{
  This is used for classic estimates, for robust estimates see \code{\link{loss.robchis}}
  This is implemented for internal use, might not be called directly by user.
}

\seealso{
\code{\link{loss.robchis}}, \code{\link{nlr}}
}
\examples{
## The function is currently defined as
"loss.chis"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ heteroscedastic variance }
\keyword{ nonlinear regression function model}
\keyword{heteroscedastic function model}
\keyword{chi-square pseudo log likelihood}
