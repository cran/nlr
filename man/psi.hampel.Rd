\name{psi.hampel}
\alias{psi.hampel}
\alias{rho.hampel}
\title{
  hampel redescending function
}
\description{
  \itemize{
    \item {psi.hampel:}  {compute hampel psi function.}
    \item{rho.hampel:} {compute hampel rho function.}
  }
}
\usage{
psi.hampel(u, a = 1.5, b = 3.5, c = 8)
rho.hampel(u, a = 1.5, b = 3.5, c = 8)
}
\arguments{
  \item{u}{
    vector of values.
}
  \item{a}{
    tuning constant, for wich less than that the quadratic function compute.
}
  \item{b}{
    tuning constant, for wich less than that the linear function compute.
}
  \item{c}{
    tuning constant, for wich less than that the constant function compute.
}
}
\details{
  psi and rho function used to in MM-estimate. The tuning constants together with k0,k1 used to give high breakdown point. This is function in MASS library, will be developed in nlr latter on.
}
\value{
  vector of psi and who function values.
}
\references{
Stromberg, A. J. (1993). Computation of High Breakdown Nonlinear Regression Parameters, Journal of American Statistical Association 88(421): 237-244.
Stromberg, J., and Ruppert, D. (1992). Breakdown in Nonlinear Regression, Journal of American Statistical Association 87: 991-997.

}
\author{
Stromberg, A. J. (1993)
}
\note{
This function is created at begining for simplicity but in feature will not used by \code{nlr} any more.
}
\seealso{
    \code{\link{rho.hampel}}
  }
\examples{
## The function is currently defined as
"psi.hampel"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Loss functions}
\keyword{Robust loss functions}
\keyword{Hampel function}
