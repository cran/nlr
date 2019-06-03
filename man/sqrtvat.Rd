\name{sqrtvat}
\alias{sqrtvat}

\title{
Compute square root of vairiance attribute.
}
\description{
Transform variance to standard deviation with all its gradient and hessian.
}
\usage{
sqrtvat(varcomp)
}

\arguments{
  \item{varcomp}{
  Is (n*1) vector of some variance, transform to \eqn{\sqrt (vc)}, with attributes attr(vc,"gradient"), \eqn{n \times p} gradient. And attr(vc,"hessian"), \eqn{n \times p \ times{p}} hessian.
}
}
\details{
For computation purpose to transform variance function values to standard deviation function value is used.
}
\value{
Standard deviation is equal the square root of variance, with Gradient equal to:
\deqn{Gradient (sdev) = \frac{1}{2} \sqrt{Var} \times Gradient (Var)} 
and hessian is equal
\deqn{
  hessian(sdev) = \frac{1}{2} \sqrt{vc} \times hesian(vc) - (\frac{1}{4} \sigma ^ 3) grad(vc)^T \%m3d\%  grad(vc)
  }
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
Is used for when standard deviation of a heterogeneous variance function model is needed.

This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{nlr}}
}
\examples{
## The function is currently defined as
"sqrtvat"
}
\keyword{Heterogeneous Variance}

