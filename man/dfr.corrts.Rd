\name{dfr.corrts}
\alias{dfr.corrts}
\title{
Derivative free Two Stage estimate
}
\description{
  Derivative free two stage estimate for nonlinear regreession model with autocorrelated error.
}
\usage{
dfr.corrts(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 0.001, minlanda = 1/2^10, 
maxiter = 25 * length(start)), correlation = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    nl.form object of the nonlinear function model. See \code{\link{nl.form}} object.
  }
  \item{data}{
    list of data with the response and predictor as name of variable.
  }
  \item{start}{
    list of starting value parameter, name of parameters must be represented as names of variable in the list.
  }
  \item{control}{
    nlr.control object, include tolerance, maxiter,... see \code{\link{nlr.control}}.
  }
  \item{correlation}{
      correlation structure, at the moment parameter of AR(p) process.
}
  \item{\dots}{
    any argument pass to \code{formula}
}
}
\details{
  In first stage nonlinear regression parameter estimate and in second stage autocorrelation structure estimate and finally the generalized least square estimates the function model parameters.
  
  In this function all stages compute by derivative free methods, which minimization methods uses Nelder-Mead method.
}
\value{
  \item{fited}{\code{\link{nl.fitt.gn}} object generated by \code{\link{nlsnm}} function.}
  \item{tm}{ fitted time series model for residuals.}

}
\references{
Riazoshams, H., Midi, H., Sharipov, O. S.H, (2010). The Performance of Robust Two Stage Estimator in Nonlinear Regression with autocorrelated Error, Communications in Statistics - Simulation and Computation, 39: 1251-1268.
}
\author{
  Hossein Riazoshams, Jul 2009.
  Email: \email{riazihosein@gmail.com}
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  This function currently run with AR process. The robust estimate is don by \code{\link{nl.robcorrts}} function. This function will be called from \code{nlr} by providing \code{correlation} as correlation structure and \code{derivfree}. It is under developement and for internal use, user might not call it directly, it is more efficient to call from \code{nlr} function with mentioned arguments.
}
\seealso{
\code{\link{nl.robcorrts}}, \code{\link{nlsqr.gn}}, \code{\link{nl.fitt.gn}}, \code{\link{nlr.control}}, \code{\link{nlsnm}}
}
\examples{
# The direct call of nlr call dfr.corrts.
p1<- 8.06e+10
p2<- 1e11
p3<-1970
p4=6
chstart2 <- list(p1=p1,p2=p2,p3=p3,p4=p4)
irandt<-nlr::trade.ir
dfrir<- dfr.corrts (nlrobj5[[4]],data=list(xr=irandt[,1],yr=irandt[,2]),start=chstart2,
control=nlr.control(trace=TRUE),correlation = 2)
dfrir$fited$parameters
}
\keyword{autocorrelated erro}
\keyword{two stage estimate}
\keyword{robust two stage estimate}
