\name{nl.hetroLS}
\alias{nl.hetroLS}
\title{
  CLSME estimate.
}
\description{
  Classic Least square based Multi Stage Estimate (CLSME) for heteroscedastic error case.
}
\usage{
nl.hetroLS(formula, data, start = getInitial(formula, data), 
control = nlr.control(tolerance = 0.001, minlanda = 1/2^10, 
maxiter = 25 * length(start)), varmodel, tau = getInitial(varmodel, vdata), ...)
}
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
  Hossein Riazoshams, May 2016, ongoing book. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  Heteroscedastic variance can have several cases, this function assume variance is parameteric function of predictor (\eqn{H(t;\tau)}). If data does not include the predictor variable of \code{varmodel} (t), the predicted of function model \eqn{f(x;\hat \theta)} will replace for (t), otherwise user have to defin (t) or (x) as predictor variable of (H).
}
\seealso{
  \code{\link{fittmethod}}, \code{\link{nl.form}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}
}
\examples{
chkft<- nl.hetroLS(formula=nlrobj1[[14]], data=list(xr=Weights$Date,yr=Weights$Weight),
start=list(p1=2200,p2=38,p3=.11), tau=list(sg=.09,landa=2),varmodel=nlrobjvarmdls1[[1]])
chkft$parameters
}
\keyword{hetroscedastic error}
