\name{optim.NM}
\alias{optim.NM}
\alias{smptry}

\title{
  NM optimization
}
\description{
  Nelder-Mead derivative free optimization. Since it is derivative free have slow convergence.
}
\usage{
optim.NM(objfnc, data, start = getInitial(objfnc, data), delta = NULL,deltar=.1, 
control=nlr.control(tolerance=1e-8, 
maxiter=250 * length(start)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{objfnc}{
    any objective function for minimizing, it must contains accept formula, data and start as argument, extra argument can be passed by (...). The output of objfnc must be a list contains: $value(attr,gradient,hessian), $angmat (angular matrix),$angvec (angular vector) to check convergence. Usually it might have \code{\link{nl.form}} object as entry.
}
  \item{data}{
    list of the data, that might have predictor and response variables with names.
}
  \item{deltar}{
    =0.1. Ratio of delta value. "nlr" compute the increament interval of golden section by start +- delta, the delta by default computed by deltar*start. The user can give direct value for deltar or value of delta to be unequally ratio movement.}
  \item{start}{
    list of initial values with names as parameters.
}
  \item{delta}{vector with same size of parameters named as parameter names, is increament for each parameter at the begining. Each parameter will be moved by delta.}
  \item{control}{
    nlr.control options to control the optimization iterations.
}
  \item{\dots}{
    any external parameters passe to \code{objfnc}.
}
}
\details{
  Optimize objective function \code{objfnc} with respect to parameters \code{start}. The mothod is derivative free using Nelder-Mead method.
  
    The function \code{smptry} Find next minimum point in Nelder-Mead algorithm. It used for internal usage might not be called by user directly.
  }
\value{
  result is a list of:
  \item{parameters}{list of estimated parameters wit hsame names as \code{start}}
  \item{objfnc}{computed object function returned back by \code{objfnc}}
  \item{history}{history of fitt, include parameters and objective values, other level of iteration is presented for which in each iteration some more steps is done to rectify the singularity of hessian.}
}
\references{
Rizo ML 2008 Statistical Computing with R The R Series. Chapman & Hall/CRC The R Series.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  User can define his own optimization function \code{objfnc} for any purpose, but this function designed eficiently for robust estimates. It is applied for minimizing several kind of objective functions such as heteroscedastic, chi-square likelihood, robust loss, but for other general problem usage is not tested. Since Nelder-Mead is derivative free its convergence is slow, so it is suggested to use more maximum number of iteration option in \code{\link{nlr.control}}. \code{\link{optim.NLM}} is more efficient method use derivative values, but when the derivative does not exist nelder-Mead can be used.

  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user. User can use it for optimization purposes.
}
\seealso{
 \code{\link{nlr.control}}, \code{\link{nl.form}}, \code{\link{optim.NLM}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"optim.NM"
}
\keyword{optimization}
\keyword{robust estimation}
\keyword{leaast square}
