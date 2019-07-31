\name{nlmest.NM}
\alias{nlmest.NM}
\alias{smptry2}
\title{
  Nonlinear MM-estimate, Nelder-Mead.
}
\description{
  MM-estimate of a nonlinear function, Using Nelder Mead derivative free optimization method.
}
\usage{
nlmest.NM(formula, data, start = getInitial(formula, data), robfunc, 
control = nlr.control(tolerance = 1e-04, 
minlanda = 1/2^25, maxiter = 100 * length(start), robscale = T), 
vm = NULL, rm = eiginv(t(chol(vm))), delta = NULL, ...)
}
\arguments{
  \item{formula}{
    nl.form object of the nonlinear function model. See \code{\link{nl.form}} object.
  }
  \item{data}{
    list of data with the response and predictor as name of variable. In heterogeneous case if it include response variable values of heterogenous variance function it asume variance function is function of predictor \eqn{H(x_i,\tau)}, otherwise it assume is a function of predictor \eqn{H(f(x_i,\theta),\tau)}.
  }
  \item{start}{
    list of starting value parameter, name of parameters must be represented as names of variable in the list.
  }
  \item{robfunc}{
    nl.form object of robust function used for downgrading.
  }
  \item{control}{
    nlr.control object, include tolerance, maxiter,... see \code{\link{nlr.control}}.
}
  \item{vm}{
    NULL, optional covariance matrix of residuals, used for nonlinear generalized M-estimate. 
}
  \item{rm}{
    optional correlation matrix, used for nonlinear generalized M-estimate. rm is correlation matrix of vm, thus only vm is enough to be given. It can be given by user also but not necessary automatically will be calculated by argument eiginv(t(chol(vm))).
  }
  \item{delta}{
    increament of Nelder Mead method, default will be calculated 10\% of parameter values, in the case of nonconvergence it can be modified manually to acheive convergence.
  }
  \item{\dots}{
    any other argument passed to formula, robfnc, or optimization function.
  }
}
\details{
  Nelder Mead is derivative free optimization method. It is used to minimize the robust loss function using \eqn{\rho} function. This method is very slow and sugest to use with a large maximum number of iterations.
  
    The function \code{smptry2} Find next minimum point in Nelder-Mead algorithm. It used for internal usage might not be called by user directly.

}
\value{
  result is object of nl.fitt.rob (nonlinear fitt robust) for homogeneous variance, and nl.fitt.rgn for heterogeneous and autocorrelated error (nonlinear fitt robust generalized), see nl.fitt.rgn object detail.
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
  if \code{vm} is not NULL the nl.fitt.rgn include following extra slots:
  \item{vm}{covariance matrix, diagonal of variance model predicted values.}
  \item{rm}{cholesky decomposition of vm.}
  \item{gresponse}{transformed of response by rm, include gradinet and hessian attributes.}
  \item{gpredictor}{transformed of predictor by rm, include gradinet and hessian attributes.}
}
\references{
  Statistical Computing with R, Maria L. Rizzo, 2008, Chopman & Hall/CRC
}
\author{
  Maria L. Rizzo
}

\note{
  This is a slow algorithm, since "nlr" is designed for derivative based, when the gradient does not exist recomend to use this function. When the gradient exists it is strongly recomend to use derivative base methods.

  This function is called from \code{\link{nlr}}, for compatibility it is more efficient to be called by \code{\link{nlr}} than callind directly.
}
\seealso{
  \code{\link{nlmest.NLM}}, \code{\link{nl.form}} 
}
\examples{

ntpstart=list(p1=.12,p2=6,p3=1,p4=33)
ntpstarttau=list(tau1=-.66,tau2=2,tau3=.04)
datalist=list(xr=ntp$dm.k,yr=ntp$cm.k)
datalist[[nlrobjvarmdls3[[2]]$independent]]<-ntp$dm.k
fittnml <- nlmest.NM(formula=nlrobj1[[15]], data = list(xr=ntp$dm.k,yr=ntp$cm.k), start=ntpstart,
robscale = TRUE, robfunc = nl.robfuncs[["hampel"]],control=nlr.control(tolerance=1e-8,trace=TRUE))
fittnml$parameters
}
\keyword{Least Square}
\keyword{Robust MM-estimate}
\keyword{Nelder-Mead optimization}

