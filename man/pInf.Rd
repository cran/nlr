\name{pInf}
\alias{pInf}
\title{
  Parameter Inference for classic nonlinear regression.
}
\description{
    Parameter inference for classic nonliner regression. It work same as \code{parInfer} method of \code{\link{nl.fitt}}, calculate covariance matrix of parameters and their confidence interval using gradient as design matrix.
}
\usage{
pInf(object, confidence = 0.95)
}
\arguments{
  \item{object}{
    Object type \code{\link{nl.fitt}} or any other of its child objects such as \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}.
}
  \item{confidence}{
    Confidence probability.
}
}
\details{
  For computing the covariance matrix of a nonlinear regression parameter, the gradient of function with respect to parameters is consider as design matrix and linear regression formulas apply for computing covariances and confidence intervals.
}
\value{
  \item{covmat:}{Covariance matrix of nonlinear model function parameters.}
  \item{corrmat}{Correlation matrix of nonlinear model function parameters.}
  \item{parstdev:}{Standard deviation of nonlinear model function parameter. It is square root of diagonal of \code{covmat}.}
  \item{CI:}{Confidence interval for nonlinear model function parameter.}
}
\references{
    Seber, G., A. F. and Wild, C. J. (2003). Nonlinear Regression. New York: John Wiley & Sons, Inc.

  Lim, C., Sen, P. K., Peddada, S. D. (2010). Statistical inference in nonlinear regression under heteroscedasticity. Sankhya B 72:202-218.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com}
  URL \url{www.riazoshams.com/nlr}
}
\note{
  This function implemented for calling for non object purpose, for example computing covarianc matrix for Weighted M-estimate stored as \code{\link{nl.fitt.rgn}} but using classic covariance computation using gradinet, instead \code{parInfer} which use convergence properties (Lim et al. 2010)
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
## The function is currently defined as
"pInf"
}
\keyword{Covariance matrix}
\keyword{Confidence interval}
\keyword{Weighted M-estimates}

