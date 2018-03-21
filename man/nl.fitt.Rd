\name{nl.fitt-class}
\Rdversion{1.1}
\docType{class}
\alias{nl.fitt}
\alias{nl.fitt-class}
\alias{$,nl.fitt-method}
\alias{[<-,nl.fitt,ANY,ANY,ANY-method}
\alias{acf,nl.fitt-method}
\alias{atypicals.deleted}
\alias{atypicals.deleted,nl.fitt-method}
\alias{atypicals}
\alias{atypicals,nl.fitt-method}
\alias{hat,nl.fitt-method}
\alias{JacobianLeverage}
\alias{JacobianLeverage,nl.fitt-method}
\alias{parInfer}
\alias{parInfer,nl.fitt-method}
\alias{plot,nl.fitt,missing-method}
\alias{predict,nl.fitt-method}
\alias{predictionI}
\alias{predictionI,nl.fitt-method}
\alias{recalc}
\alias{recalc,nl.fitt-method}
\alias{residuals,nl.fitt-method}
\alias{plot,nl.fitt-method}
\title{Class \code{"nl.fitt"}}
\description{
  Object for a fitted nonlinear regression model. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nl.fitt", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{parameters}:}{Object of class \code{"list"}, estimate of nonlinear model \eqn{\theta}.}
    \item{\code{scale}:}{Object of class \code{"numericorNULL"}, standard deviation scale estimate \eqn{\sigma}.}
    \item{\code{correlation}:}{Object of class \code{"numericorNULL"}, correlation structure of error.}
    \item{\code{form}:}{Object of class \code{"nl.form"} of nonlinear model.}
    \item{\code{response}:}{Object of class \code{"vectororMatrix"} response, left side of formula.}
    \item{\code{predictor}:}{Object of class \code{"vectororMatrix"}, estimated predictor \eqn{\eta(\hat{\theta)}}.}
    \item{\code{curvature}:}{Object of class \code{"listorNULL"} of PE and IE curvatures.}
    
    \item{\code{history}:}{Object of class \code{"matrixororNULL"} convergence computations in iteration procedures, include parameters, objective function and other parameters depends on the method.}
    
    \item{\code{method}:}{Object of class \code{"fittmethodorNULL"} method of iteration used, contains main method, functions and sub methods. See \code{\link{fittmethod}}.}
    
    \item{\code{data}:}{Object of class \code{"list"} data used in computation, including response and predictor variables.}
    
    \item{\code{sourcefnc}:}{Object of class \code{"callorNULL"} source function called for fitt.}
    
    \item{\code{Fault}:}{Object of class \code{"Fault"} of error or warnings if happened.}
    
    \item{\code{others}:}{Object of class \code{"listorNULL"} of other computations, as an example the object of outlier detection measures will be saved in this slot later on.}
  }
}
\section{Extends}{
Class \code{"\link{nl.fittorNULL}"}, directly.
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "nl.fitt")}: access the slot values. Usage: \code{objectname$slotname}}
    
    \item{atypicals}{\code{signature(nlfited = "nl.fitt")}: detect atypical points by calculating outlier detection measures. Usage: \code{atypicals(nlfited)}}
    
    \item{hat}{\code{signature(x = "nl.fitt")}: generalized Hat matrix from linear regression to nonlinear regression using gradient. Usage: \code{hat(x="nl.fitt")}}
    
    \item{JacobianLeverage}{\code{signature(nlfited = "nl.fitt")}: Jacobian-Leverage for nonlinear regression. Usage \code{JacobianLeverage(nlfited = "nl.fitt")}}

    \item{parInfer}{\code{signature(object = "nl.fitt")}: parameter inference function, calculate covariance matrix of parameters and their confidence interval using gradient as design matrix. Usage: \code{parInfer(object,confidence = .95)}
    }
    \item{plot}{\code{signature(x = "nl.fitt", y = "missing",control=nlr.control(),...)}: generic function extended to \code{nl.fitt} object. Plot the object. Usage. \code{plot(x,y="missing",control=nlr.control(history=F,length.out=NULL,singlePlot=F))}. If history is TRUE the convergence of fitt will be ploted. length.out is length of incremented p[redictor to acheive smooter curve. singlePlot=F plot the model and residuals in two collumn.}
    
    \item{predict}{\code{signature(object = "nl.fitt")}: generic function, predict nonlinear function model at estimated parameter values. Usage: \code{predict(object,...)}, dots argument can include the \code{newdata} which might be new  list of new values for predictor variables, if not given the original data that used for fitt will be used to calculate prediction values.}

    \item{predictionI}{\code{signature(nlfited = "nl.fitt.gn")}: prediction interval. Usage: \code{predictionI(nlfited,confidence=.95,data=NULL)}, data is new data that will be predicting the values for them.}
    
    \item{recalc}{\code{signature(object = "nl.fitt")}: generic function, recalculate the object with new arguments given in dots argument.}
    
    \item{residuals}{\code{signature(object = "nl.fitt")}: residuals of fitt. Usage: \code{residuals((object,...))}, dots argument can include \code{data} list of predictor and response variables, if \code{data} is not given the residuals will calculate for original data used in estimation.}
    \item{acf}{ The function acf computes (and by default plots) estimates of the autocovariance or autocorrelation function of residuals. For argument details see stats  \code{\link[=stats]{acf}} general function.}
	 }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  All information of a nonlinear fited model are saved in \code{nl.fitt}, thus it can be large variable of informations. The generalized form \code{\link{nl.fitt.gn}} and robust forms \code{\link{nl.fitt.rob}} and \code{\link{nl.fitt.rgn}} of a fitt is children of this object. Typically it used to save Least-Square estimation method. But it is extensively used to save other fited objects such as hetroscedastic variance parameter fitts.
}
\seealso{
  \code{\link{nl.fitt.gn}}, \code{\link{fittmethod}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}, \code{\link{Fault}}, \code{\link{nl.form}}, \code{\link[=stats]{acf}}
}
\examples{
showClass("nl.fitt")
}
\keyword{classes}
\keyword{robust fitt}
\keyword{least square}
\keyword{heteroscedastic variance}
