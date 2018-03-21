\name{optim.NLM}
\alias{optim.NLM}
\title{
  NLM optimization.
}
\description{
  Modified Newton-Levenberg-Marquardt optimization. It is derivative based optimization method, designed to be robust against sigularity problem due to outliers.
}
\usage{
optim.NLM(objfnc, data, start = getInitial(objfnc, data), 
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
  Optimize objective function \code{objfnc} with respect to parameters \code{start}. The mothod is gradient base combines Newton, Stepest descend and levenberg-Marquardt. 
}
\value{
  result is a list of:
  \item{parameters}{list of estimated parameters wit hsame names as \code{start}}
  \item{objfnc}{computed object function returned back by \code{objfnc}}
  \item{history}{history of fitt, include parameters and objective values, other level of iteration is presented for which in each iteration some more steps is done to rectify the singularity of hessian.}
}
\references{
Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
    Seber, G., A. F. and Wild, C. J. (2003). Nonlinear Regression. New York: John Wiley & Sons, Inc.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  User can define his own optimization function \code{objfnc} for any purpose, but this function designed eficiently for robust estimates. It is applied for minimizing several kind of objective functions such as heteroscedastic chi-square likelihood, robust loss, but for other general problem usage is not tested.
  
  This function call by \code{\link{nlr}}, for compatibility it is better to call from \code{\link{nlr}} rather than directly by user. User can use it for optimization purposes.
}
\seealso{
  \code{\link{nl.form}}
}
\examples{
## The function is currently defined as
"optim.NLM"
}
\keyword{optimization}
\keyword{robust estimation}
\keyword{leaast square}