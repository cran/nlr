\name{atyps}
\alias{atyps}
\title{
 Find atypical points.}
\description{
  This is extension of \code{\link{nlout}} for which add the jacobian leverage to the output. It is designed for feature development.
}
\usage{
  atyps(nlfited)
}

\arguments{
  \item{nlfited}{
    An \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}} obejct after a nonlinear model fitted before.
  }
}
\details{
  Compute outlier detection measures by calling \code{\link{nlout}} function and the ncompute jacobian leverages then dispatches the output.
}
\value{
  List of jacobian leverage and list of outlier detection measures derived by  \code{\link{nlout}}.
  \itemize{
    \item{jlev}{\code{jlev} jacobian leverages.}
    \item{nlout}{list of output from   \code{\link{nlout}} function.}
    
  }
}
\references{
  Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  Used mostly for internal purposes.
}
\seealso{
  \code{\link{nlsqr}}
}
\examples{
#Example
#d<-list(xr=Weights$Date, yr=Weights$Weight)
# wmodel <- nlr(nlrobj1[[2]],data=d,control=nlr.control(method = "OLS",trace=TRUE))
 #a=nlout(wmodel)
#atyps(wmodel)
}
\keyword{QR Decomposition}
\keyword{Full Rank QR decomposition}
