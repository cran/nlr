\name{dfr.hetroLS}
\alias{dfr.hetroLS}
\title{
  Derivative free CLSME.
}
\description{
  Derivative free Classic Least square based Multi Stage Estimate (CLSME) for heteroscedastic error case.
}
\usage{
dfr.hetroLS(formula, data, start = getInitial(formula, data), control = nlr.control(
tolerance = 1e-04, minlanda = 1/2^10,maxiter = 25 * length(start)), varmodel, 
tau = getInitial(varmodel, vdata), ...)
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
  Least square based estimate for nonlinear regression with hetroscedastic error when variance is a general function of unkown parameters.
}
\value{
  generalized  fitt object \code{\link{nl.fitt.gn}}. The \code{hetro} slot include parameter estimate and other information of fitt for heteroscedastic variance model.
  \item{(parameters }{nonlinear regression parameter estimate of \eqn{\theta}.}
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
  \item{vm}{covariance matrix, diagonal of variance model predicted values.}
  \item{rm}{cholesky decomposition of vm.}
  \item{gresponse}{transformed of response by rm, include gradinet and hessian attributes.}
  \item{gpredictor}{transformed of predictor by rm, include gradinet and hessian attributes.}
  \item{hetro}{
    \code{\link{nl.fitt}} object of fited variance odel:
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
      }
  }
}
\references{
  Riazoshams, H. (2012), Robustifying the Least Squares estimate of parameters of variance model function in nonlinear regression with heteroscedastic variance, Poster Presentation, Royal Statistical Society Conference (RSS) 2012, Telford, UK.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{
  Heteroscedastic variance can have several cases, this function assume variance is parameteric function of predictor (\eqn{H(t;\tau)}). If data does not include the predictor variable of \code{varmodel} (t), the predicted of function model \eqn{f(x;\hat \theta)} will replace for (t), otherwise user have to defin (t) or (x) as predictor variable of (H).
  
    \code{dfr.hetroLS} is derivative free it is slow convergence, while \code{\link{nl.hetroLS}} is derivative based estimate is effectively fast method. Since it is slow algorithm it is recomneded to use larger values for maximum number of iterations in \code{\link{nlr.control}} options.
}
\seealso{
  \code{\link{fittmethod}}, \code{\link{nl.form}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.hetroLS}}, \code{\link{nlr.control}}
}
\examples{
ntpstart=list(p1=.12,p2=6,p3=1,p4=33)
ntpstarttau=list(tau1=-.66,tau2=2,tau3=.04)
datalist=list(xr=ntp$dm.k,yr=ntp$cm.k)
htls<- dfr.hetroLS(formula=nlrobj1[[15]], data=datalist, start= ntpstart,tau=ntpstarttau,
varmodel=nlrobjvarmdls3[[2]],control=nlr.control(tolerance=1e-8))
htls$parameters
}
\keyword{hetroscedastic error}
\keyword{Nelder-Mead optimization}

