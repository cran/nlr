\name{nl.robhetro}
\alias{nl.robhetro}
\title{
  Robust Multi Stage Estimate.
}
\description{
  (RME) for nonlinear regression with heteroscedastic variance, when the variance of error is general parameteric function of unkown parameters. Robust form of CME (See \code{\link{nl.hetro}}).
}
\usage{
nl.robhetro(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 1e-05, minlanda = 1/2^10, 
maxiter = 25 * length(start)), robfunc, varmodel, tau = NULL, ...)
}
\arguments{
  \item{formula}{
    \code{\link{nl.form}} object of the nonlinear function model.
}
  \item{data}{
    list of data include responce and predictor.
}
  \item{start}{
    list of parameter values of nonlinear model function (\eqn{\theta}.
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
  In stage 1 the nonlinear model parameter estimates by robust MM-estimate, Stage 2 compute robust sample variance of data, Stage 3 estimate the parameter of variance function model by maximizing the robustified form of chi-square pseudo-likelihood function. Stage 4 estimate the final value of function model parameter by generalized robust MM-estimate.
}
\value{
  \code{\link{nl.fitt.rgn}} for heterogeneous and autocorrelated error (nonlinear fitt robust generalized) will return.
  \item{parameters }{nonlinear regression parameter estimate of \eqn{\theta}.}
  \item{correlation}{of fited model.}
  \item{form}{\code{\link{nl.form}} object of called nonlinear regression model.}
  \item{response}{computed response.}
  \item{predictor}{computed (right side of formula) at estimated parameter with gradient and hessian attributes.}
  \item{curvature}{list of curvatures, see \code{\link{curvature}} function.}
  \item{history}{matrix of convergence history, collumns include: convergence index, parameters, minimized objective function, convergence criterion values, or other values. These values will be used in \code{plot} function in ploting history.}
  \item{method}{\code{\link{fittmethod}} object of method used for fitt.}
  \item{data}{list of called data.}
  \item{sourcefnc}{Object of class \code{"callorNULL"} source function called for fitt.}
  \item{Fault}{\code{\link{Fault}} object of error, if no error Fault number = 0 will return back.}
  \item{htheta}{robust loss value including gradient and hessain attributes.}
  \item{rho}{computed robust rho function, including gradient and hessain attributes.}
  \item{ri}{estimated residuals, including gradient and hessain attributes.}
  \item{curvrob}{curvature}
  \item{robform}{\code{nl.form} object of robust loss rho function.}
  \item{vm}{covariance matrix, diagonal of variance model predicted values.}
  \item{rm}{cholesky decomposition of vm.}
  \item{gresponse}{transformed of response by rm, include gradinet and hessian attributes.}
  \item{gpredictor}{transformed of predictor by rm, include gradinet and hessian attributes.}
  \item{hetro}{
    \code{\link{nl.fitt.rob}} object of fited variance odel:
      \itemize{
        \item{parameters}{estimate of variance parameter \eqn{\tau}}
        \item{form}{\code{nl.form} object of called \code{varmodel}.}
        \item{predictor}{variance model computed at estimated parameter, \eqn{H(x;\hat\tau)}}
        \item{response}{sample variance computed used as response variable.}
        \item{history}{matrix of convergence history, collumns include: convergence index, parameters, minimized objective function, convergence criterion values, or other values.}
        \item{method}{\code{\link{fittmethod}} object of method used for fitt.}
        \item{data}{response (\eqn{z_i}) and predictor \code{t} variable values, used to computing the variance model.}
        \item{sourcefnc}{Object of class \code{"callorNULL"} source function called for fitt.}
        \item{Fault}{\code{\link{Fault}} object of error, if no error Fault number = 0 will return back.}
  			\item{htheta}{robust loss value including gradient and hessain attributes, for variance model. In fact is loglikelihood values.}
        \item{rho}{computed robust rho function, including gradient and hessain attributes.}
      }
  }
  \item{others}{$refvar reference variance. variance of zi's.}
}
\references{
  Riazoshams, H,. 2010. Outlier detection and robust estimation methods for nonlinear regression having autocorrelated and heteroscedastic errors. PhD thesis disertation, University Putra Malaysia.
}
\author{
Hossein Riazoshams
}
\note{
  Heteroscedastic variance can have several cases, this function assume variance is parameteric function of predictor (\eqn{H(t;\tau)}). If data does not include the predictor variable of \code{varmodel} (t), the predicted of function model \eqn{f(x;\hat \theta)} will replace for (t), otherwise user have to defin (t) or (x) as predictor variable of (H).
}
\seealso{
  \code{\link{fittmethod}}, \code{\link{nl.form}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
# ntp data fitt
# tolerance is set as 1e-3 for testing purposes
# is not accurate enough, user can increase it.
ntpstart=list(p1=.12,p2=6,p3=1,p4=33)
ntpstarttau=list(tau1=-.66,tau2=2,tau3=.04)
datalist=list(xr=ntp$dm.k,yr=ntp$cm.k)
datalist[[nlrobjvarmdls3[[2]]$independent]]<-ntp$dm.k
aa1 <- nl.robhetro(formula=nlrobj1[[16]],data=datalist,
start=ntpstart,robfunc=nl.robfuncs[["hampel"]],
tau=ntpstarttau,varmodel=nlrobjvarmdls3[[2]],robscale=T,method="NM",
control=nlr.control(tolerance=1e-4))
aa1$parameters
}
\keyword{robust statistics}
\keyword{heteroscedastic error}
\keyword{MM-estimators}