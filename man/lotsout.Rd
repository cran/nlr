\name{lotsout}
\alias{lotsout}
\docType{data}
\title{
Artificially Contaminated Data from Logistic model. 
}
\description{
lotsout is artificially Contaminated Data from Logistic model, using the computed parameter values from the estimates of fited chicken growth \code{\link{Weights}} data.
}
\usage{data("lotsout")}
\format{
  The format is:
 chr "lotsout"
  \itemize{
  \item lotsout[,1]: predictors repeated two times.
  \item lotsout[,2]: response simulated from the parameter estimates and logistic model.
 }
}
\details{
The predictor is similar to predictor of chicken growth \code{\link{Weights}} data but only repeated two times. The response is simulated from logistic model using the parameters from estimated fit for chicken growth data. This selection is used to mimic a real data example to have more natural behaviour.
}
\source{
 Robust Nonlinear Regression, with Application Using R, Riazoshams et al, Wiley Inc.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Jan 2010. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\examples{
lotsout
}
\keyword{datasets}
