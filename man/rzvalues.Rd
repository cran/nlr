\name{rzvalues}
\alias{rzvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Robust sample variance
}
\description{
  Compute Robust Sample variance for dependent variable. For repeated response the sample Mead Absolute Deviance (MAD)  of predictor compute over a cross section of predictor, for non repeated data sample variance of response computed by difference of residuals of consequetive values.
}
\usage{
rzvalues(res, ni, xo)
}
\arguments{
  \item{res}{
    The data to compute variance for. In most of application it is residuals of fitt.
}
  \item{ni}{
    vector of bumber of repeated data. It can be output from \code{\link{nonrepl}} function.
}
  \item{xo}{
    Position of the repeated data in original vector.
}
}
\details{
  Typically it is used to compute the robust variance of residuals output from \code{\link{nonrepl}} function.
}
\value{
  vector of robust sample variance.
}
\references{
  Riazoshams, H. (2012), Robustifying the Least Squares estimate of parameters of variance model function in nonlinear regression with heteroscedastic variance, Poster Presentation, Royal Statistical Society Conference (RSS) 2012, Telford, UK.
   
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
  }
\author{
    Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com}
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  This is robustified form of sample variance function \code{\link{zvalues}}.
}
\seealso{
   \code{\link{zvalues}}
}
\examples{
## The function is currently defined as
"rzvalues"
}
\keyword{sample variance}
\keyword{Robust variance}
