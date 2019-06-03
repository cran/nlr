\name{transfquad}
\alias{transfquad}
\alias{transfquadvec}
\title{
  Variance to standard deviation transform.
}
\description{
  This function compute the standard deviation of an output from a parameteric nonlinear variance function object, usualy stored as \code{\link{nl.form}}. The gradient and hessian attributes also will be transformed.
  \code{transfquadvec} transform a vector similarly.
}
\usage{
transfquad(varcomp)
}
\arguments{
  \item{varcomp}{
     Is vector of variances, with (n by p) \code{"gradient"} and (three dimentional n by p by p) \code{"hessian"} attributes.
}
}
\details{
  The standard deviation is simply square root of variance. The gradinet is transformed of square root of gradient. Hessian is transformed of second derivative of square root function as well.
}
\value{
  Vector of transformed standard deviation, including (n by p) \code{"gradient"} and (three dimentional n by p by p) \code{"hessian"} attributes.
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
  In \code{nlr} variance function is stored in \code{\link{nl.form}} for computations.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{nl.form}}
}

\examples{
## The function is currently defined as
"transfquad"
}
\keyword{ variance }
\keyword{ standard deviation }
\keyword{gradient}
\keyword{hessian}
