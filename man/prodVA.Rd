\name{prodVA}
\alias{prodVA}
\title{
  Product of three dimentional array in vector.  
}
\description{
  Array(n*p*p)*Vector(n*1) is (p*p) matrix equal the sum of vector multiple first dimention of array.
}
\usage{
prodVA(ary, vector)
}
\arguments{
  \item{ary}{
    (n*p*P) numeric array.
}
  \item{vector}{
    (p*p) matrix.
}
}
\details{
  It can be used to multiple a (n*p*p) Hessian into (n*1) vector of residuals.
}
\value{
  (p*p) matrix.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  It is used in optimization and loss functions, may not explicitly called by user. May be used in defining loss functions by user.
}
\seealso{
  \code{\link{\%3d*m\%}}, \code{\link{\%m3d\%}}
}
\examples{
## The function is currently defined as
"prodVA"
a=array(1:8,c(2,2,2))
b=c(2,3)
d=prodVA(a,b)
}
\keyword{Hessian}
\keyword{Gradient}% __ONLY ONE__ keyword per line
