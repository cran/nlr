\name{nl.robhetroWM}
\alias{nl.robhetroWM}
\title{
  Weighted M-estimate.
}
\description{
  Weighted M-estimate is robustified form of MLE, for nonlinear regression with heteroscedastic error, when variance is parameteric fucntion form. Both nonliner regression model parameter and variance function parameters compute simultanously by minizing the robustified objective function form.
}
\usage{
nl.robhetroWM(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 1e-04, minlanda = 1/2^10, 
maxiter = 50 * length(start), derivfree = T), robfunc, varmodel, tau = varmodel$par, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
      \code{\link{nl.form}} object of the nonlinear function model.
}
  \item{data}{
      list of data include responce and predictor.
}
  \item{start}{
      list of parameter values of nonlinear model function (\eqn{\theta}. in \eqn{f(x,\theta)}).
}
  \item{control}{
      list of \code{\link{nlr.control}} for controling convergence criterions.
}
  \item{robfunc}{
    nl.form object of robust function used for downgrading.
}
  \item{varmodel}{
      \code{nl.fomr} object of variance function model for heteroscedastic variance.
}
  \item{tau}{
      list of initial values for variance model function \code{varmodel} argument.
}
  \item{\dots}{
      extra arguments to nonlinear regression model, heteroscedastic variance function, robust loss function or its tuning constants.
}
}
\details{
  For minimizing the objective function simultaneously for theta and tau, derivative free method Nelder-Mead is used.
}
\value{
  return object \code{\link{nl.fitt.rgn}} for nonlienar regression with heterogeneous error.  
  \item{parameters }{nonlinear regression parameter estimate of \eqn{\theta}.}
  \item{correlation}{of fited model.}
  \item{form}{\code{\link{nl.form}} object of called nonlinear regression model.}
  \item{response}{computed response.}
  \item{predictor}{computed (right side of formula) at estimated parameter with gradient and hessian attributes.}
  \item{history}{matrix of convergence history, collumns include: convergence index, parameters, minimized objective function, convergence criterion values, or other values. These values will be used in \code{plot} function in ploting history.}
  \item{method}{\code{\link{fittmethod}} object of method used for fitt.}
  \item{data}{list of called data.}
  \item{sourcefnc}{Object of class \code{"callorNULL"} source function called for fitt.}
  \item{Fault}{\code{\link{Fault}} object of error, if no error Fault number = 0 will return back.}
  \item{htheta}{robust loss likelihood value including gradient and hessain attributes.}
  \item{rho}{computed robust rho function, including gradient and hessain attributes.}
  \item{ri}{estimated residuals, including gradient and hessain attributes.}
  \item{robform}{\code{nl.form} object of robust loss rho function.}
  \item{vm}{covariance matrix, diagonal of variance model predicted values.}
  \item{rm}{cholesky decomposition of vm.}
  \item{hetro}{
    \code{\link{nl.fitt.rob}} object of fited variance odel:
      \itemize{
        \item{parameters}{estimate of variance parameter \eqn{\tau}}
        \item{form}{\code{nl.form} object of called \code{varmodel}.}
        \item{predictor}{variance model computed at estimated parameter, \eqn{H(x;\hat\tau)}}
        \item{response}{sample variance computed used as response variable.}
      }
  }
  \item{others}{$refvar reference variance. variance of zi's.}
}
\references{
Lim, C., Sen, P. K., Peddada, S. D. (2010). Statistical inference in nonlinear regression under heteroscedasticity. Sankhya B 72:202-218.
}
\author{
  Lim, C., Sen, P. K., Peddada, S. D. 
}

\note{
  Heteroscedastic variance can have several cases, this function assume variance is parameteric function of predictor (\eqn{H(t;\tau)}). If data does not include the predictor variable of \code{varmodel} (t), the predicted of function model \eqn{f(x;\hat \theta)} will replace for (t), otherwise user have to defin (t) or (x) as predictor variable of (H).
}
\seealso{
 \code{\link{fittmethod}}, \code{\link{nl.form}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
ntpstart=list(p1=.12,p2=6,p3=1,p4=33)
ntpstarttau=list(tau1=-.66,tau2=2,tau3=.04)
datalist=list(xr=ntp$dm.k,yr=ntp$cm.k)
datalist[[nlrobjvarmdls3[[2]]$independent]]<-ntp$dm.k
# ntp data fitt
# tolerance is set as 1e-3 for testing purposes
# is not accurate enough, user can increase it.
bb1 <- nl.robhetroWM(formula=nlrobj1[[15]],data=datalist,
start=ntpstart,robfunc=nl.robfuncs[["least square"]],
tau=ntpstarttau,varmodel=nlrobjvarmdls3[[2]],control=nlr.control(tolerance=1e-3,maxiter=1500))
bb1$parameters
#---------------- hampel -----------------
aa1 <- nl.robhetroWM(formula=nlrobj1[[15]],data=datalist,start=ntpstart,
robfunc=nl.robfuncs[["hampel"]],derivfree=T,
tau=ntpstarttau,varmodel=nlrobjvarmdls3[[2]],
control=nlr.control(tolerance=1e-3,maxiter=1500))#,delta=c(0.2,1,1,160,.2,1,.03))
aa1$parameters
}
\keyword{robust statistics}
\keyword{heteroscedastic error}
\keyword{Weighted M-estimators}