\name{nlout.JL}
\alias{nlout.JL}
\title{
  Nonlinear outlier detection.
}
\description{
  Detecting outlier for nonlinear regression, is based on mixing statsitics measures and robust estimates through their covariance matrices (hat matrix). The covariance matrix in nonlinear is based on the gradient of nonlinear regression model, but it based on linear approximation of the model, instead Jacobian Leverage is used in this function.
}
\usage{
nlout.JL(nlfited)
}
\arguments{
  \item{nlfited}{
  Object of type \code{\link{nl.fitt}} or \code{\link{nl.fitt.gn}} for classic estimators, \code{\link{nl.fitt.rob}} or \code{\link{nl.fitt.rgn}} for robust estimators.
}
}
\details{
  The outlier detection measutred used in this function are studentized residuals and Cook Distance. They are mixture of estimators and Jacobians. They are successful for detecting outlier only if combine with robust fits, eventhough the function can work with classic fits but it is not recomended. 
  Riazoshams and Midi (2014) 
}
\value{
  Result is list of \code{\link{nl.robmeas}} objects for each statistics. 
  \item{"jl.vmat"}{
      Jacobian-leverage matrix.
      }
  \item{"jl.studres"}{ \code{\link{nl.robmeas}} object of Jacobian Leverage Studentised Residuals.}
  \item{"jl.cook"}{ \code{\link{nl.robmeas}} object of Jacobian Leverage Elliptic Norm (Cook Dist).}
  \item{"jl.hadi"}{\code{\link{nl.robmeas}} object of Jacobian Leverage Hadi potential.}
  \item{"jl.delstud"}{ \code{\link{nl.robmeas}} object of Jacobian Leverage Deletion Studentized.}
  \item{"jl.dffits"}{ \code{\link{nl.robmeas}} object of Jacobian Leverage DFFITS.}
  \item{"jl.atk"}{ \code{\link{nl.robmeas}} object ofJacobian Leverage Atkinson Distance.}
}
\references{
    Riazoshams H, Habshah M and Adam MB 2009 On the outlier detection in nonlinear regression. 3(12), 243-250.
  
  Riazoshams H and Midi H 2014 Robust Leverage and outlier detection measures in nonlienar regression, 2014 (Unpublished manuscript).
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  This function return back all resutls and statistics based on Jacobian leverage, but Riazoshams (2014) showed studentized residuals when combine with robust estimators can detect outliers correctly. Thus to identify outlier correctly first estimate the parameters bu robust options of \code{\link{nlr}} function then call \code{nlout}, finally look at the list values \code{"jl.delstud"} from the result list. 
  The \code{plot} and other methods of  \code{\link{nl.robmeas}} display the results visually. 
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{nl.robmeas}}, \code{\link{nlr}}, \code{\link{nlout}}
}
\examples{
d<-list(xr=Weights$Date, yr=Weights$Weight)
wmodel <- nlr(nlrobj1[[2]],data=d,control=nlr.control(method = "OLS",trace=TRUE))
a=nlout.JL(wmodel)
plot(a[[2]])
}
\keyword{ Mahalanobis distance }
\keyword{ Cook Distance}
\keyword{ Jacobian Leverage }
\keyword{ Studentized Residuals}
