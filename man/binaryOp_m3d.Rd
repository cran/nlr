\name{\%m3d\%}
\alias{\%m3d\%}
\title{
  Matrix 3 dimentional product.
}
\description{
  Binary Operator, for sample size (n), vector with size p of gradient stored in (n*p) matrix, transpose of gradient multiple by gradient for (n) points can be obtained by (gradient \%m3d\% gradient). Result save in (n*p*p) three dimentional array.
}
\usage{
  mat1 \%m3d\%  mat2
}
\arguments{
  \item{mat1}{
    n*p matrix, in gradient example n is sample size, p is number of parameters, mat1[i,] is (p) vector of gradient of ith sample.
}
  \item{mat2}{
    p*n matrix.
}
}
\details{
    It is Used to implement gradient product for n sample 
}
\value{
  \deqn{Array_{n*p*p} =\nabla_{\theta}^T f_i(\theta) \nabla_{\theta} f_i(\theta),i=1,\dots,n}
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\seealso{
  \code{\link{\%3d*m\%}}, \code{\link{prodVA}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"\%m3d\%"
m1 <- matrix(c(1:4),nrow=2)
m2 <- matrix(c(1:4),nrow=2)
d=m1 \%m3d\% m2
}
\keyword{ Hessian}
\keyword{ Gradient}% __ONLY ONE__ keyword per line

