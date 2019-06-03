\name{mscale}
\alias{mscale}
\title{
  Scale M-estimate
}
\description{
  Scale M-estimator with 50\% breakdown
}
\usage{
mscale(u)
}
\arguments{
  \item{u}{
    vector of values, in nonlinear regression residuals are used.  
  }
}
\details{
  Robust M-estimate of scale using rodust loss function rho, to acheive 50\% breakdown. Tunning constant k1 should be stored in \code{\link{nl.form}} object of the robust function.
  
  This is called from \code{\link{nlmest.NLM}} function to derive MM.estimate.
}
\value{
Minimized Sum of rho function, by reweighting.
}
\references{
Yohai (1987) Annals, Stromberg (1993) JASA. GKS  2 June 99
}
\author{
Stromberg (1993)
}
\note{
    General form is developed in \code{nlr} package to work for general robust function.
  }
\seealso{
  \code{\link{nl.form}}, \code{\link{nlmest.NLM}}
}
\examples{
x=c(2,3,10)
mscale(x)
}
\keyword{Scale}
\keyword{Scale M-estimate}
\keyword{MM-estimate}
