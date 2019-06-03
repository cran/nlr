\name{transform}
\alias{transform}

\title{
Transform by R matrix.
}
\description{
Transform \eqn{R \times value} and its gradient and hessian.
}
\usage{
transform(value,rm)
}

\arguments{
  \item{value}{
A vector to be transformed by R-matrix. It must have gradient and hessian as attributes.
}
  \item{rm}{
  The R-matrix, from QR decomposition.
}
}
\details{
  It is used to transfomr a nonlinear regression model by R matrix from QR decomposition of variance covariance matrix.
}
\value{
 Is a transformed vector \eqn{R \times V}, and hessian and gradient as its attributes. Gradient is equal:
 \deqn{gradient= rm \times g}
 and hessian equal the three dimentional product:
 \deqn{hessian= h \%3d*m\% rm}.
}
\references{
  Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
It is used in generalized model. In both heteroscedastic and autocorrelated model will be applied.

  This function call by \code{\link{nlr}}, it might not be called directly by user.
}

\seealso{
\code{\link{nlmest.NLM}}, \code{\link{nlmest.NLMf}}, \code{\link{nlmest.WF}}
}
\examples{
## The function is currently defined as
"transform"
}

\keyword{Generalized Least Square}
\keyword{Generalized Nonlinear Regression model}
