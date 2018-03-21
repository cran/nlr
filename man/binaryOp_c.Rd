\name{\%c\%}
\alias{\%c\%}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  cross product binary operator.
}
\description{
  Binary operator, given matrices x and y as arguments, return a matrix cross-product.
}
\usage{
  x \%c\%  y
}
\arguments{
  \item{x}{
    (m*n) matrix.
}
  \item{y}{
    (m*q) matrix.
}
}
\details{
    result is (n*m) matrix = t(x) \%*\% y, sum(x[,i]*y[,j])
}
\value{
  A double or complex matrix, with appropriate dimnames taken from x and y.}
\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  This is binary operator of \code{crossprod} adjusted in \code{nlr} for more streamlined and readable. Mostly used in optimization functions, might not be called by user explicitly.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
(1:4) \%c\% (1:4) ## result is sum(1^2+2^2+3^2+4^2)
 m1 <- matrix(c(1:4),nrow=2)
 m2 <- matrix(c(1:8),nrow=2)
 m1 \%c\% m2
}
\keyword{Matrix multiplication}

