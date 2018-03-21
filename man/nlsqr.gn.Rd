\name{nlsqr.gn}
\alias{nlsqr.gn}
\title{
  Generalized Least Square estimate.
}
\description{
  Generalized Least Square estimate of a nonlinear function, Using QR-decomposition of Gradient matrix.
}
\usage{
nlsqr.gn(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 0.001,
minlanda = 1/2^10, maxiter = 25 * length(start)), vm, rm = eiginv(t(chol(vm))))
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
    Covariance matrix of residuals, used for nonlinear generalized M-estimate. 
}
  \item{rm}{
    optional correlation matrix, used for nonlinear generalized M-estimate. rm is correlation matrix of vm, thus only vm is enough to be given. It can be given by user also but not necessary automatically will be calculated by argument eiginv(t(chol(vm))).
  }
}
\details{
  It is used to minimize the square loss function, using QR-decomposition of gradient matrix, thus the nonlinear function model \code{formula} must return back Gradient. \code{nlsqr.gn} work with a general variance covariance matrix, such as heteroscedastic or weights in variance, and partially autocorrelated with any general format.
}
\value{
  result is object of \code{\link{nl.fitt.gn}} for generalized fitt when covariance or correlation matrix is given.
  \item{parameters}{nonlinear regression parameter estimate of \eqn{\theta}.}
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
}
\references{
  Seber, G., A. F. and Wild, C. J. (2003). Nonlinear Regression. New York: John Wiley & Sons, Inc.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com}
  URL \url{www.riazoshams.com/nlr}
}
\note{
  This is a generalized form of \code{\link{nlsqr}} function. It is fast algorithm based on gradient. If gradient does not exist one can use \code{\link{nlsnm}} function.

  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user.
}
\seealso{
  \code{\link{nl.form}}, \code{\link{nlsnm}}, \code{\link{nlr.control}}, \code{\link{nl.fitt}}, \code{\link{curvature}}, \code{\link{Fault}},\code{\link{nlsqr}}
}
\examples{
## The function is currently defined as
"nlsqr.gn"
}
\keyword{Least Square}
\keyword{Gradient Matrix}
\keyword{optimization}
\keyword{Generalized Least Square}