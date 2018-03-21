\name{zvalues}
\alias{zvalues}
\title{
  Sample Variance
}
\description{
  Compute Sample variance for dependent variable. For repeated response the sample variance of predictor compute over a cross section of predictor, for non repeated data sample variance of response computed by difference of residuals of consequetive values.
}
\usage{
zvalues(res, ni, xo)
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
  Typically it is used to compute the variance of residuals output from \code{\link{nonrepl}} function.
}
\value{
  vector of classic variance.
}
\references{
Bunke, O., Droge, B., Polzehl, J. Splus tools for model selection in nonlinear regression (1998) Computational Statistics, 13 (2), pp. 257-281.
}
\author{
Bunke, O., Droge, B., Polzehl.
}
\note{
  The robustified form of this function is \code{\link{rzvalues}}.
  This function call by \code{\link{nlr}}, it might not be called directly by user.
}

\seealso{
 \code{\link{rzvalues}}
}
\examples{
## The function is currently defined as
"zvalues"
}
\keyword{sample variance}
\keyword{Robust variance}
