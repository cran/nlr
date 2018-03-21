\name{nlmest.RWT}
\alias{nlmest.RWT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonlinear MM-estimate using reweighting method.
}
\description{
Compute MM-estimate using reweighting method developed by Stromberg.
}
\usage{
nlmest.RWT(formula, data, start = getInitial(formula, data), robfunc, 
control = nlr.control(tolerance = 0.001, minlanda = 1/2^25, 
maxiter = 25 * length(start), trace = F), vm = NULL, rm = eiginv(t(chol(vm))), ...)
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
  \item{vm}{optional covariance matrix of residuals, used for nonlinear generalized M-estimate.
}
  \item{rm}{optional correlation matrix, used for nonlinear generalized M-estimate. rm is correlation matrix of vm, thus only vm is enough to be given. It can be given by user also but not necessary automatically will be calculated by argument eiginv(t(chol(vm))).
}
  \item{\dots}{any other argument passed to formula, robfnc, or optimization function.
  }
  \item{control}{\code{nlr.control} option variables.}
}
\details{
Compute MM-estimate using reweighting method developed by Stromberg.
}
\value{
  result is object of nl.fitt.rob (nonlinear fitt robust) for homogeneous variance, and nl.fitt.rgn for heterogeneous (not developed yet) and autocorrelated error (nonlinear fitt robust generalized), see nl.fitt.rgn object detail.
}
\references{
Stromberg, A. J. (1993). Computation of High Breakdown Nonlinear Regression Parameters, Journal of American Statistical Association 88(421): 237-244.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
It is similar to \code{nlrob} function in \code{robustbase} package, it is not completly operational in \code{nlr}.
  
  This function is called from \code{\link{nlr}}, for compatibility it is more efficient to be called by \code{\link{nlr}} than callind directly.
}
\seealso{
\code{\link{nlmest.NLM}}
}
\examples{
## The function is currently defined as
"nlmest.RWT"
}
\keyword{MM-estimate}
\keyword{Iteratively reweighting}
