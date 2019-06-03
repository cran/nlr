\name{individ}
\alias{individ}
\title{
  Split individuals
}
\description{
  Hessian of a function is symetric matrix, \code{individ} function constract lower diagonal parts and store them in collumns of a matrix. Thus the only unique hessians will be given.
}
\usage{
individ(hessian)
}
\arguments{
  \item{hessian}{
    Three dimentional \eqn{n \times p \times p} array of hessian.
  }
}
\details{
  It is used in \code{\link{curvature}} function to compute the curvatures.
  }
\value{
  \eqn{n \times (p*p+1)} matrix of hessian values, constracted from the lower triangular of hessian.
}
\references{
Bates, D.M., and Watts, D. G. (1980). Relative curvature measures of nonlinearity, J. R. statistic. Ser. B 42: 1-25.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}}
\note{
  It is created for internal use in \code{\link{curvature}} function, it ight not be called by user explicitly.
}

\seealso{
\code{\link{curvature}}
}
\examples{
## The function is currently defined as
"individ"
}
\keyword{curvature}
\keyword{gradient}
\keyword{hessian}

