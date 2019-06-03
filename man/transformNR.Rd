\name{transformNR}
\alias{transformNR}
\alias{transforminv}
\title{
  Transform nonlinear regression model
}
\description{
  Transform both side a nonlinear regression model by cholesky decomposition of covariance matrix of errors. The new transformed model has constant variance and uncorrelated errors.
  \itemize{
    \item{\code{transform}:}{ Transform using \code{rm} matrix}
    \item{\code{transforminv}:}{ Inverse transform using inverse of \code{rm} matrix}
  }
}
\usage{
transformNR(value, rm)
}
\arguments{
  \item{value}{
    Numeric value with (n by p) \code{"gradient"} and (three dimentional n by p by p) \code{"hessian"} attributes. Usually it is right  side (predictors) or left side (response) of a nonlinear regression model.
}
  \item{rm}{
    Matrix for transforming. In nonlinear regression if covaraince matrix of error represent by \eqn{\sigma ^2 V}, \code{rm} is cholesky decompostion of \code{V}.
}
}
\details{
  If errors of a nonlinear regression is Autocorrelated or heteroscedastic, the model can be transformed to uncorrelated and homoscedastic using cholesky decomposition.
  \itemize{
    \item \code{transformNR}: Multiply by \code{rm} matrix, create standardized model.
    \item \code{transforminv}: Multiply by inverse of \code{rm} matrix, create heteroscedastic or autocorrelated model.
  }
}
\value{
    Transformed values including (n by p) \code{"gradient"} and (three dimentional n by p by p) \code{"hessian"} attributes.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.

Seber, G., A. F. and Wild, C. J. (2003). Nonlinear Regression. New York: John Wiley & Sons, Inc.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  For transfering a nonlinear regression model both side should be transformed. In \code{nlr} package generalized nonlinear model use this transformation.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{transforminv}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"transformNR"
}
\keyword{heteroscedastic variance function}
\keyword{Autocorrelated error}
\keyword{Cholesky decomposition}

