\name{prodAV}
\alias{prodAV}
\alias{prodVAV}
\title{
Compute product of Array into a Vector.
}
\description{
function: prodAV product array to matrix
}
\usage{
prodAV(ary, vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ary}{
  is (p*p*n) array.
}
  \item{vector}{
  is a vector with dimension (p)
}
}
\details{
  function: prodAV product array to matri.
  compute: 'array * vector', which is (n,p) vector.

}
\value{
'array * vector', which is (n,p) vector.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
Three dimentional array multiple matrix have several possibilities, see refeerence. 
\code{prodVAV} is another case of vector multiple array multiple vector. Is for feature purposes.
}
\seealso{
  \code{\link{\%3d*m\%}}, \code{\link{\%m3d\%}} \code{\link{prodVA}}
}
\examples{
## The function is currently defined as
"prodAV"
}
\keyword{Three dimentional matrix product}
\keyword{Three dimentional array}

