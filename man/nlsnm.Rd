\name{nlsnm}
\alias{nlsnm}
\title{
  Least Square estimate.}
\description{
  Least Square estimate of a nonlinear function, Using Nelder Mead derivative free optimization method.
}
\usage{
  nlsnm(formula, data, start = getInitial(formula, data), delta = NULL, 
  control = nlr.control(tolerance = 1e-04, 
  maxiter = 100 * length(start)), vm = NULL, rm = NULL, ...)
}
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
  Nelder Mead is derivative free optimization method. It is used to minimize the square loss function. This method is very slow and sugest to use with a large maximum number of iterations.
}
\value{
  result is object of \code{\link{nl.fitt}} (nonlinear fitt robust) for homogeneous variance, and \code{\link{nl.fitt.gn}} for generalized fitt when covariance or correlation matrix is given.
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
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{nlmest.NLM}}, \code{\link{nl.form}}, \code{\link{nlsqr}}
}
\examples{
## The function is currently defined as
"nlsnm"
}
\keyword{Least Square}
\keyword{Nelder-Mead optimization}
\keyword{optimization}