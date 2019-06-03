\name{jaclev}
\alias{jaclev}
\title{Jacobian Leverage for nonlinear regression.
}
\description{ Compute the Jacobian Leverage, generalized for nonlinear case.
}
\usage{
jaclev(gradient, hessian, rsd)
}
\arguments{
  \item{gradient}{\eqn{n \times p} gradient of nonlinear function.
}
  \item{hessian}{three simentional \eqn{n \times p \times p} of hessian of nonlinear regression function.
}
  \item{rsd}{\eqn{n \times 1} residual vector.
}
}
\details{
  Jacobian leverage, generalized form of hat matrix for nonlinear regression.
}
\value{
  \eqn{n \times n} matrix of jacobian leverages.
}
\references{Laurent. R. T. ST., and Cook. R. D. (1992). Leverage and Superleverage in Nonlinear Regression, Journal of the American Statistical Association 87(420): 985-990.
}
\author{
Laurent. R. T. ST., and Cook.
}
\note{
  Jacobian leverage for nonlinear regression is direct definition of perturbing response, thus it is free from the problems due to linear approximation of nonlinear function.
}

\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}
  }
\examples{
## The function is currently defined as
"jaclev"
}
\keyword{ Jacobian Leverage}
