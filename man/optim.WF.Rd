\name{optim.WF}
\alias{optim.WF}
\alias{CubInterp}
\alias{zoom}
\title{
  WF optimization
}
\description{
  Optimization using Wolf conditions.
}
\usage{
optim.WF(objfnc, data, start = getInitial(objfnc, data), 
control = nlr.control(tolerance = 0.001, minlanda = 1/2^10, 
maxiter = 25 * length(start)), ...)
}
\arguments{
  \item{objfnc}{
    any objective function for minimizing, it must contains accept formula, data and start as argument, extra argument can be passed by (...). The output of objfnc must be a list contains: $value(attr,gradient,hessian), $angmat (angular matrix),$angvec (angular vector) to check convergence. Usually it might have \code{\link{nl.form}} object as entry.
}
  \item{data}{
    list of the data, that might have predictor and response variables with names.
}
  \item{start}{
    list of initial values with names as parameters.
}
  \item{control}{
    nlr.control options to control the optimization iterations.
}
  \item{\dots}{
    any external parameters passe to \code{objfnc}.
}

}
\details{
  Optimize objective function \code{objfnc} with respect to parameters \code{start}. The mothod is gradient base using Wolf condition for rectifying the negative definit hessian problems. 
  
  The following function are called from \code{optim.WF}.
  
  \code{CubInrep}: Cubic interpolation, 
}
\value{
  result is a list of:
  \item{parameters}{list of estimated parameters wit hsame names as \code{start}}
  \item{objfnc}{computed object function returned back by \code{objfnc}}
  \item{history}{history of fitt, include parameters and objective values, other level of iteration is presented for which in each iteration some more steps is done to rectify the singularity of hessian.}
}
\references{
Nocedal J and Wright SJ 2006 Numerical optimization.. New York, NY.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  User can define his own optimization function \code{objfnc} for any purpose, but this function designed eficiently for robust estimates. It is applied for minimizing several kind of objective functions such as heteroscedastic chi-square likelihood, robust loss, but for other general problem usage is not tested.
  
  Together with \code{\link{optim.NLM}} are used in \code{nlr} package to acheive optimization result in failure of one another method.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user. User can use it for optimization purposes.
}
\seealso{
  \code{\link{nl.form}}, \code{\link{optim.NLM}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"optim.WF"
}
\keyword{optimization}
\keyword{robust estimation}
\keyword{leaast square}
\keyword{Wolf conditions}
