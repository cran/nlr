\name{nlr}
\alias{nlr}
\title{Non-Linear Robust fitt.
}
\description{This is generic function fits a nonlinear mixed-effects model using robust methods described in Riazoshams et al (2015), allowed to be correlated and/or have unequal variances.
}
\usage{
nlr(formula, data = parent.frame(), start = getInitial(formula, data), 
control = nlr.control(minlanda=1 / 2 ^ 10, 
maxiter=25 * length(start)), weights = NULL, robustobj = NULL, 
robustform = c("hampel", "huber", "bisquare",
"andrew", "halph huber", "hampel 2", "least square"), 
varianceform = NULL, tau = NULL, correlation = NULL, 
covariance = NULL, ...)
}
\arguments{
  \item{formula}{ 
  can be a nl.form object of the nonlinear function model, or a model formula with the response on the left of a ~ operator and an expression involving parameters and covariates on the right. See \code{\link{nl.form}} object.
}
  \item{data}{an optional data frame or list of data with the response and predictor as name of variable. In heterogeneous case if it include response variable values of heterogenous variance function it asume variance function is function of predictor \eqn{H(x_i,\tau)}, otherwise it assume is a function of predictor \eqn{H(f(x_i,\theta),\tau)}.
}
  \item{start}{list of starting value parameter, name of parameters must be represented as names of variable in the list.
}
  \item{control}{nlr.control, include controles for nlr, see \code{\link{nlr.control}}.
}
  \item{weights}{User optional matrix of variance covariance matrix of error, a general weight that can be used by user, but if correlation, covariace or varianceform argument given, then will be ignored.
}
  \item{robustobj}{=NULL, optional \code{nl.form} object of robust loss function defined by user. User can define his/her own function, to be nl.form object, Or use defined robloss function in robustform.
}
  \item{robustform}{nl.form object of robust function used for downgrading. functions("hampel","huber","bisquare","andrew","halph huber","hampel 2","least square").
}
  \item{varianceform}{NULL, \code{nl.form} object if given heterogeneous variance will be fited. 
}
  \item{tau}{NULL, list or data.frame of initial values for heterogeneous variance function parameter. the stored value in vardnc object of nl.form will be stored.
}
  \item{correlation}{autocorrelated error, form of corStruct but not in nlme, is extra arguments to any of "nl.forms".
}
  \item{covariance}{optional covariance matrix of errors. If given Generalized estimates will be calculated.
}
  \item{\dots}{Any extra arguments to any function sources, such as nonlinear regression model, heteroscedastic variance function, robust loss function or optimization object function.
}
}
\details{nlr is non-linear robust inference. It is mostly gradient and hessian based. The classic estimates also is avilable.
}
\value{
  Output depends on the method called as:
  \itemize{
    \item {Nonlinear Least Square Estimate(NLLS):} {\code{\link{nl.fitt}} object.}
    \item {Generalized NLLS:} {\code{\link{nl.fitt.gn}} object.}
    \item{Robust Estimate:} {\code{\link{nl.fitt.rob}} object.}
    \item {Generalized Robust Estimate:} {\code{\link{nl.fitt.rgn}} object.}
  }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{
  Thic is the global function of \code{nlr} package. You can fitt models both with classic and robust methods. Besides homogeneous variance of error, heteroscedastic and autocorrelated error can be fited by this function.
  
  \code{nlr} is optimized for derivative based computation, but deribative free methods are provided for nonlinear regressions that derivative does not exist. In contrast derivative free method in this package are slow.
  
  This package requires \code{tseries} package to be installed, under ubuntu if facing with error that can not install "curl" ot "libcurl" package, run these command in terminal:
  sudo apt-get install libcurl4-openssl-dev    

deb http://security.ubuntu.com/ubuntu/ precise-security restricted main multiverse universe
deb http://us.archive.ubuntu.com/ubuntu/ precise-updates restricted main multiverse universe

sudo apt-get update
sudo apt-get install libcurl4-gnutls-dev

after all we can install tseries in R.
}
\seealso{
  \code{\link{nl.form}}, \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
## Chicken fitt without hetroscedasticity
 d<-list(xr=Weights$Date, yr=Weights$Weight)
 fitt.chicken1 <- nlr(nlrobj1[[2]],data=d,control=nlr.control(method = "OLS",trace=TRUE))
fitt.chicken1$parameters
plot(fitt.chicken1)
## Chicken fitt with hetroscedasticity
# RME by default
fitt.chicken2<-nlr(formula=nlrobj1
                   [[14]],data=d,start=list(p1=2300,p2=42,p3=.11), 
                   robustform ="hampel",
                   tau=list
                   (sg=.09,landa=1),varianceform=nlrobjvarmdls1[[1]],
                   control=nlr.control(tolerance=1e-3))
fitt.chicken2$parameters
fitt.chicken2$hetro$parameters
plot(fitt.chicken2)
# autocorrelated case
xr = trade.ir[, 1]
yr = trade.ir[, 2]
a1 <-  nlr( nlrobj5[[18]],data=list(xr = xr, yr = yr),
correlation=list(StructName="corAR1"))
a1$parameters
plot(a1)
}
\keyword{Nonlinear regression}
\keyword{Robust statistics}
\keyword{Robust Nonlinear regression}
