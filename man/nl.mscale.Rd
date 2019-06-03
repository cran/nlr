\name{nl.mscale}
\alias{nl.mscale}
\title{
Scale M-estimator with 50\% breakdown.
}
\description{
Compute High Breakdown point M-estimate of scale \eqn{\sigma}.
}
\usage{
nl.mscale(u, robfunc, ...)
}

\arguments{
  \item{u}{
Residuals \eqn{r_i=y_i-f(x_i;\theta)}.
}
  \item{robfunc}{
Robust \eqn{\rho_2} function used to compute M-estimate of scale/
}
  \item{\dots}{
  Ane other parameter passed to \eqn{\rho} function, and others.
}
}
\details{
This estimate is used in MM-estimate procedure of location parameter \eqn{\theta}.
}
\value{
  Single Numeric value of \eqn{\sigma} estimate.
}
\references{
  Yohai (1987) Annals, Stromberg (1993) JASA.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  Its value is used inside MM-estimate procedure of parameter \eqn{\theta}. Currently work with Hampel \eqn{\rho} function, not all of \eqn{\rho} defined functions.
  It is called by \code{\link{nlmest.NLM}} in scale estimation steps, might not be called directly by user.
}
\seealso{
\code{\link{nlmest.NLM}}
}
\examples{
## The function is currently defined as
"nl.mscale"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MM-estimate}
\keyword{Scale Estimate}
\keyword{Location parameter estimate}
