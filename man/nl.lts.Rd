\name{nl.lts}
\alias{nl.lts}

\title{
Compute (LTS) Least Trimmed Square Estimate.
}
\description{
LTS is minimizing trimmed sum of squares. 
}
\usage{
nl.lts(formula, data, start, h = NULL, control = nlr.control())
}
\arguments{
  \item{formula}{
    \code{\link{nl.form}} object of the nonlinear model function.
}
  \item{data}{
  List of data used in predictor and response.
}
  \item{start}{
  List of parameter starting value.
}
  \item{h}{
    Percentage of trimming the residuals, if omited the default 25\% will be used.
}
  \item{control}{
  \code{\link{nlr.control}} options, will be submited to least square \code{\link{nlsqr}} function.
}
}
\details{
  LTS trimme \code{h} percent of residuals first then compute the least square estimate, and final parameter estimate is the one minimize the sum of squares of errors.
}
\value{
  \code{\link{nl.fitt}} object of fitted model parameter.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
The result data returnd in fitted object is trimmed data.
}

\seealso{
 \code{\link{nl.form}}, \code{\link{nl.fitt}}
}
\examples{
data=list(xr=Weights$Date,yr=Weights$Weight)
fit<- nl.lts(nlrobj1[[14]],data=data,start=list(p1=1000,p2=42,p3=.11))
fit$parameters
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Least Square Estimae}
\keyword{Least Trimmed Square}
