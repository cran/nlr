\name{nl.form-class}
\Rdversion{1.1}
\docType{class}
\alias{nl.form}
\alias{nl.form-class}
\alias{$,nl.form-method}
\alias{$<-,nl.form-method}
\alias{[<-,nl.form,ANY,ANY,ANY-method}
\alias{all.vars,nl.form-method}
\alias{eval,nl.form-method}
\alias{eval.nl.form}
\alias{evald,nl.form-method}
\alias{getInitial,nl.form-method}
\alias{selfStart,nl.form-method}

\title{Class \code{"nl.form"}}
\description{It is a class of nonlinear regression model function, can embed expression formula and function. It is a more general object such that the heterogeneous variance model and robust loss function in \code{nlr} is saved as \code{nl.form} object.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form:

\code{new("nl.form",formula,fnc,formtype,p, inv=NULL,name=name,par,arguments=list(...),dependent, independent,origin,selfStart) }, 

or constructor:

\code{nl.form(form,p=NULL,inv=NULL,name,par=NULL,dependent=NULL,independent=NULL,origin=NULL,selfStart=NULL,...)}
}
\section{Slots}{
  \describe{
    \item{\code{formula}:}{Object of class \code{"call or NULL"} it can be 1- a two sided formula with response  (or a function of response) in left of ~ and nonlinear function model is a function of predictors and parameters, or 2- a one sided formula with ~nonlinear model in right, is again a function of predictors and unknown parameters. If Null then the nonlinear model is a R function stored in \code{fnc} slot as bellow. Each of right side or left side formula can return "gradient" and "hessian" as attributed value.}
    \item{\code{fnc}:}{Object of class \code{"function or NULL"} is nonlinear model stored as R function.}
    \item{\code{dependent}:}{Object of class \code{"characterorNULL"} of predictor variable, null value means the \code{formula} slot is one sided.}
    \item{\code{independent}:}{Object of class \code{"characterorNULL"} caharacter vector name of predictor variables, which can be more than one predictor.}
    \item{\code{formtype}:}{Object of class \code{"character"} character type of \code{nl.form}. Do not insert this slot, it will be set automatically by creater, it use internaly for evaluation.}
    \item{\code{p}:}{Object of class \code{"numericorNULL"} number of parameters.}
    \item{\code{inv}:}{Object of class \code{"callorNULL"} if nonlinear model is function of one predictor the inverse function define in this slot. Still not functioning, designed for feature extention, so it can be ignored at the moment, but it is better to define.}
    \item{\code{name}:}{Object of class \code{"character"} a character name for the nonlinear function model.}
    \item{\code{par}:}{Object of class \code{"list"} of parameters, assigned value to parameters will be used in worst case that initial values can not be computed.}
    \item{\code{arguments}:}{Object of class \code{"list"} list of extra arguments use in \code{formula} or \code{fnc} slot.}
    \item{\code{origin}:}{Object of class \code{"callorNULL"} is original one sided or two sided expression of nonlinear model withought gradient and hessian.}
    \item{\code{selfStart}:}{Object of class \code{"functionorNULL"} selfstart function defined for initial values guiss. Is same as nls or nlme::nlme functions.}
  }
}
\section{Extends}{
Class \code{"\link{nl.formorNULL}"}, directly.
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "nl.form")}: return slots.}
    \item{all.vars}{\code{signature(expr = "nl.form")}: Return a character vector containing all the names which occur in \code{"formula"} or or \code{"fnc"}.}
    \item{eval}{\code{signature(expr = "nl.form")}: evaluate \code{nl.form} object in the environment include parameters and predictor variables.}
    \item{evald}{\code{signature(expr = "nl.form")}: same functionality as eval but extended for compatibility. \code{eval.nl.form} has same functionality created for compatibility purpose.}
    \item{getInitial}{\code{signature(object = "nl.form")}: get initial value from selfstart or par slot.}
    \item{selfStart}{\code{signature(model = "nl.form")}: self computing initial value, if not given \code{"getInitial"} method return initial values from environment or "par" slot.}
    \item{eval.nl.form}{\code{signature(expr = "eval.nl.form")}: \code{eval.nl.form} has same functionality as \code{eval}, created for compatibility purpose.}
	 }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons. 
}
\author{
  Hossein Riazoshams, 2013. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  this object typicaly implemented to store a nonlinear regression model function informations. But extensively used in \code{nlr} package to save heteroscedastic variances and robust loss functions.
}
\seealso{
  \code{\link{nl.robfuncs}}
}
\examples{
# define hampel robust loss as a function
hampel <- nl.form(
	  form = function(t,
	                  a = 1.345,
	                  k0 = 3.73677,
	                  k1 = 4,
	                  maxrho5 = 1.345,
	                  ...) {
	    U <- abs(t)
	    Ugrta <- (U > abs(a))
	    .rho <- .grad <- .hess <- .weight <- NULL
	    .rho[Ugrta] <- 2. * abs(a) * U[Ugrta] - a * a
	    .rho[!Ugrta] <- t[!Ugrta] ^ 2
	    .grad[Ugrta] <- 2. * abs(a) * sign(t[Ugrta])
	    .grad[!Ugrta] <- 2. * t[!Ugrta]
	    .hess[Ugrta] <- 0.
	    .hess[!Ugrta] <- 2.
	    .weight[Ugrta] <- 2. * abs(a) / U[Ugrta]
	    .weight[!Ugrta] <- 2.
	    attr(.rho, "gradient") <- .grad
	    attr(.rho, "hessian") <- .hess
	    attr(.rho, "weight") <- .weight
	    return(.rho)
	  },
	  name = "huber",
	  independent = "t",
	  a = 1.345,
	  k0 = 3.73677,
	  k1 = 4,
	  maxrho5 = 1.345
	)
}
\keyword{classes}
\keyword{nl.form}
