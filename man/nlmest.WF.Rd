\name{nlmest.WF}
\alias{nlmest.WF}
\alias{zoom2}
\title{
Nonlinear MM-estimate using wolf conditions.
}
\description{
  Parameters estimates by robust MM-estimate by minimizing the sum of robust rho function, Choosing a Step-Length using Wolfe Conditions.
}
\usage{
nlmest.WF(formula, data, start = getInitial(formula, data), robfunc, 
control = nlr.control(tolerance = 1e-04, 
maxiter = 25 * length(start), robscale = T), vm = NULL, rm = eiginv(t(chol(vm))), ...)
# zoom2 is sub function
#zoom2(a1,a2,p1,p2,pd1,pd2,ht,phi0,phid0,theta1,
#delta1,sigma,objfnc,data,start,control,...)
}
\arguments{
  \item{formula}{nl.form object of the nonlinear function model. See \code{\link{nl.form}} object.
}
  \item{data}{list of data with the response and predictor as name of variable. In heterogeneous case if it include response variable values of heterogenous variance function it asume variance function is function of predictor \eqn{H(x_i,\tau)}, otherwise it assume is a function of predictor \eqn{H(f(x_i,\theta),\tau)}.
}
  \item{start}{list of starting value parameter, name of parameters must be represented as names of variable in the list.
}
  \item{robfunc}{nl.form object of robust function used for downgrading.
}
  \item{control}{nlr.control object, include tolerance, maxiter,... see \code{\link{nlr.control}}.
}
  \item{vm}{optional covariance matrix of residuals, used for nonlinear generalized M-estimate.
}
  \item{rm}{optional correlation matrix, used for nonlinear generalized M-estimate. rm is correlation matrix of vm, thus only vm is enough to be given. It can be given by user also but not necessary automatically will be calculated by argument eiginv(t(chol(vm))).
}
  \item{\dots}{any other argument passed to formula, robfnc, or optimization function.
}
}
\details{
  This function Choosing a Step-Length using Wolfe Conditions in direct search optimization. Sum of robust loss function \eqn{\rho} is minimized in order to get the robust MM-estimates. This method is very fast and used for when the gradient of colaborating functions exists. 
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
Numerical Optimization, Jorge Nocedal Stephen J. Wright, Springer 2006.
}
\author{
Jorge Nocedal Stephen J. Wright
}
\note{
The \code{\link{nlmest.NLM}} function is more developed form of this function, but still there are some situations that current function might be used. \code{nlr} package for estimating robust MM, try to acheive the optimization convergence using \code{\link{nlmest.NLM}} or \code{\link{nlmest.WF}} or \code{\link{nlmest.LM}}, either of them that does not attain convergence other method might be used.
}
\seealso{
\code{\link{nlmest.LM}},  \code{\link{nlmest.NLM}},  \code{\link{nlmest.LM}}, \code{\link{nlmest.WF}}, \code{\link{nl.form}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{nlr.control}}
}
\examples{
## The function is currently defined as
"nlmest.WF"
}
\keyword{Robust statistics}
\keyword{MM-estimate}
\keyword{Optimization}
\keyword{Levenberg-Marquardt}
\keyword{Newton optimization}
\keyword{direct search}
\keyword{wolfe conditions}


