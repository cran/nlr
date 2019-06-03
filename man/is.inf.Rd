\name{is.inf}
\alias{is.inf}
\alias{is.missing}
\title{
Check IEEE Arithmetic Values}
\description{
Returns a logical vector or matrix describing the type of numeric elements present. This distinguishes between infinite values, NaN's, missing values and ordinary numbers. 
}
\usage{
is.inf(x)
}
\arguments{
  \item{x}{
    numeric vector of values. Check performs for all values. 
  }
}

\details{
These functions are created for compatibility with SPLUS.
}
\value{
Returns an object similar to the input which is filled logical values. Values will be false for vectors that are not of mode "numeric". 

is.infinite is TRUE for values of x that are either plus or minus infinity. 

is.inf is an abbreviation for is.infinite 
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, 2013. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
    In \code{nlr} package this function use in \code{\link{eval}} method of \code{\link{nl.form}} to identify the computed response and predictor is right.
  }
\seealso{
  \code{\link{nl.form}}
}
\examples{
aaray=c(0,3/0)
is.inf(aaray)
is.missing(aaray)
}
\keyword{infinite}

