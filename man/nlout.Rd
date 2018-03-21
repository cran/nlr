\name{nlout}
\alias{nlout}
\title{
  Nonlinear outlier detection.
}
\description{
  Detecting outlier for nonlinear regression, is based on mixing robust estimates and statsitics measures.
}
\usage{
nlout(nlfited)
}
\arguments{
  \item{nlfited}{
    Object of type \code{\link{nl.fitt}} or \code{\link{nl.fitt.gn}} for classic estimators, \code{\link{nl.fitt.rob}} or \code{\link{nl.fitt.rgn}} for robust estimators.
}
}
\details{
  The outlier detection measutred used in this function are studentized residuals and Cook Distance. They are mixture of estimators and Jacobians. They are successful for detecting outlier only if combine with robust fits, eventhough the function can work with classic fits but it is not recomended. 
}
\value{
  Result is list of \code{\link{nl.robmeas}} objects for each statistics. 
  \item{"vmat"}{
      variance covariance matrix of parameters
      \deqn{\sigma^2 (\nabla f(\theta)'\nabla f(\theta))^{-1})}
      }
  \item{"d.yhat"}{ predicted values after rremoving a point \eqn{\hat y_{(-i)}}}
  \item{"studres"}{ \code{\link{nl.robmeas}} object Studentized residuals.}
  \item{"cook"}{ \code{\link{nl.robmeas}} object od Elliptic Norm (Cook Dist)}
  \item{"mahd.v"}{\code{\link{nl.robmeas}} object of Regression Mahalanobis Distance.}
  \item{"mahd.dt"}{ \code{\link{nl.robmeas}} object of Mahalanobis MVE, data.}
  \item{"mahd.xs"}{ \code{\link{nl.robmeas}} object of  Mahalanobis MVE, xs.}
	\item{"hadi"}{ \code{\link{nl.robmeas}} object of Hadi potential.}
	\item{"potmah"}{ \code{\link{nl.robmeas}} object of Potential mahalanobis.}
  \item{"delstud"}{ \code{\link{nl.robmeas}} object of Deletion Studentized.}
  \item{"dffits"}{ \code{\link{nl.robmeas}} object of DFFITS.}
	\item{"atk"}{ \code{\link{nl.robmeas}} object of Atkinson Distance.}
	\item{"mvedta"}{ \code{\link{nl.robmeas}} object of MVE data.}
  \item{"mvex"}{ \code{\link{nl.robmeas}} object of MVE x.}
	\item{"dfbetas"}{ \code{\link{nl.robmeas}} object of DFBETAS.}
}
\references{
    Riazoshams H, Habshah M and Adam MB 2009 On the outlier detection in nonlinear regression. 3(12), 243-250.
}
\author{
  Hossein Riazoshams, Dec 2008
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  This function return back all resutls and statistics but, Riazoshams (2009) showed studentized residuals and Cook distance when combine with robust estimators can detect outliers correctly. Thus to identify outlier correctly first estimate the parameters bu robust options of \code{\link{nlr}} function then call \code{nlout}, finally look at the list values \code{"studres"} and \code{"cook"} from the result list. 
  The \code{plot} and other methods of  \code{\link{nl.robmeas}} display the results visually. 
}

\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{nl.robmeas}}, \code{\link{nlr}}, \code{\link{nlout.JL}}
}
\examples{
 d<-list(xr=Weights$Date, yr=Weights$Weight)
 wmodel <- nlr(nlrobj1[[2]],data=d,control=nlr.control(method = "OLS",trace=TRUE))
 a=nlout(wmodel)
 plot(a[[3]]) # is in fact wrong inference due to hetroscedasticity
}
\keyword{ Mahalanobis distance }
\keyword{ Cook Distance}
\keyword{ Studentized Residuals}
