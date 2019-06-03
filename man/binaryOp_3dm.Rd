\name{\%3d*m\%}
\alias{\%3d*m\%}
\title{
  Product array to matrix
}
\description{
Binary operator, compute multiplication of three dimentional array to a vector, sum over first dimention of array to the vector. Technically cross section values in first dimention are equal. It use for computation purposes in optimizations.
}
\usage{
  ary \%3d*m\% vector
}
\arguments{
  \item{ary}{
    (n*p*p) numeric array.
}
  \item{vector}{
    (n*1) numeric vector.
}
}
\details{
  Three dimentional array multiple by vector, first dimention of array and vector related to sample data, multiple of them and sum is a fixed value will be returned back.
}
\value{
  result(:,i,j) = sum (ary[,i,j]*vector), but repeat in all cross section, thus the cross section values are equal.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2013. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{
  It use for computation purposes in optimization, might not be called directly by user.
}
\seealso{
\code{\link{\%m3d\%}}, \code{\link{prodVA}}
}
\examples{
  ##---- Should be DIRECTLY executable !! ----
  ##-- ==>  Define data, use random,
  ##--	or do  help(data=index)  for the standard data sets.
  ## The function is currently defined as
  "\%3d*m\%"
  a=array(1:8,c(2,2,2))
  b=3:4
  d=a \%3d*m\% b
}
\keyword{ three dimentional array }
