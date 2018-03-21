\name{nl.robfuncs}
\alias{nl.robfuncs}
\docType{data}
\title{
  Robust Loss functions provided for \code{nlr}.
}
\description{
  List of nl.form objects of 7 pre defined robust loss functions. Robust loss functions used for robust estimating parameters. They can be used in several part of the package functions.
}
\usage{nl.robfuncs}
\format{
  The format is: list
 nl.robfuncs[index]
  \itemize{
    \item \code{nl.robfuncs[1]:} huber function.
    \item \code{nl.robfuncs[2]:} hampel function.
    \item \code{nl.robfuncs[3]:} bisquare function.
    \item \code{nl.robfuncs[4]:} andrew function.
    \item \code{nl.robfuncs[5]:} halph huber function.
    \item \code{nl.robfuncs[6]:} hampel 2 function.
    \item \code{nl.robfuncs[7]:} least square (quadratic) function.
 }
}
\details{Each of loss functions include tuning parameters as extra argument to \code{\link{nl.form}}. The result is the computed loss function with attributes of "gradient", "hessian", "weights". the object are mostly stored in fnc slot which is function, it can be called directly as a function.
}
\source{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\examples{
data(nl.robfuncs)
## maybe str(nl.robfuncs) ; plot(nl.robfuncs) ...
## the object are mostly stored in fnc slot which is function, it can be called
## directly as a function.
plot(seq(-6,6,length.out=30),nl.robfuncs[[1]]$fnc(seq(-6,6,length.out=30)),type="l",
xlab="t",ylab="rho",main=nl.robfuncs[[1]]$name)
}
\keyword{datasets}
