\name{dfr.hetro}
\alias{dfr.hetro}
\title{
  Derivative free (CME)
}
\description{
  (CME) Classic multi stage estimate for nonlinear regression with heteroscedastic error, when variance is function of unkown parameters. The variance function model parameter estimate using pseudo chi-square likelihood of computed sample variance. \code{dfr.hetro} is derivative free estimate CME.
}
\usage{
dfr.hetro(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 1e-05, minlanda = 1/2^10, 
maxiter = 25 * length(start)), varmodel, tau = NULL, ...)
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
  
  \item{varmodel}{
      \code{nl.fomr} object of variance function model for heteroscedastic variance.
  }
  
  \item{tau}{
      list of initial values for variance model function \code{varmodel} argument.
  }
  \item{\dots}{
      extra arguments to nonlinear regression model, heteroscedastic variance function, robust loss function or its tuning constants, or optimization functions.
  }
}
\details{
  In stage 1 the nonlinear model parameter estimates by Classic OLS, Stage 2 compute sample variance of data, Stage 3 estimate the parameter of variance function model by maximizing the chi-square pseudo-likelihood function. Stage 4 estimate the final value of function model parameter by generalized least square. For optimization the derivative free Nelder-Mead is used.
}
\value{
  generalized  fitt object \code{\link{nl.fitt.gn}}. The \code{hetro} slot include parameter estimate and other information of fitt for heteroscedastic variance model.
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
Riazoshams, H,. 2010. Outlier detection and robust estimation methods for nonlinear regression having autocorrelated and heteroscedastic errors. PhD thesis disertation, University Putra Malaysia.

   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  Heteroscedastic variance can have several cases, this function assume variance is parameteric function of predictor (\eqn{H(t;\tau)}). If data does not include the predictor variable of \code{varmodel} (t), the predicted of function model \eqn{f(x;\hat \theta)} will replace for (t), otherwise user have to defin (t) or (x) as predictor variable of (H).
  
  \code{dfr.hetro} is derivative free it is slow convergence, while \code{\link{nl.hetro}} is derivative based estimate is effectively fast method.
}
\seealso{
  \code{\link{fittmethod}}, \code{\link{nl.form}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.hetro}}
}
\examples{
ntpstart22=list(p1=.12,p2=7,p3=1,p4=160)
ntpstarttau22=list(tau1=-1.24,tau2=2.56,tau3=.03042)
datalist=list(xr=ntp$dm.k,yr=ntp$cm.k)
datalist[[nlrobjvarmdls3[[2]]$independent]]<-ntp$dm.k
ntpfit<- dfr.robhetro(formula=nlrobj1[[16]],data=datalist,start=ntpstart22,
robfunc=nl.robfuncs[["hampel"]], tau=ntpstarttau22,
varmodel=nlrobjvarmdls3[[2]],robscale=TRUE,method="NM",control=nlr.control(tolerance=1e-4,
maxiter=150))
ntpfit$parameters
}
\keyword{hetroscedastic error}
\keyword{Nelder-Mead optimization}
