\name{nl.fitt.rgn}
\Rdversion{1.1}
\docType{class}
\alias{nl.fitt.rgn}
\alias{parInfer,nl.fitt.rgn}
\alias{nl.fitt.rgn-class}
\alias{parInfer,nl.fitt.rgn-method}
\alias{predictionI,nl.fitt.rgn-method}
\alias{residuals,nl.fitt.rgn-method}
\alias{atypicals,nl.fitt.rgn-method}

\title{Class \code{"nl.fitt.rgn"}}
\description{
  Object of generalized robust estimates of nonlinear regression model.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nl.fitt.rgn", ...)}.
}
\section{Slots}{
  \describe{
    Robust generalized slots:
    \item{\code{vm}:}{Object of class \code{"matrix"} of variance covariance matrix of error.}
    
    \item{\code{rm}:}{Object of class \code{"matrix"} of correlated error.}
    
    \item{\code{hetro}:}{Object of class \code{"nl.fittorNULL"}, include object \code{nl.fitt} of heteroscedastic error fit, or \code{NULL} for non hetroscedastic. It include parameter estimates of hetroscedastic variance \eqn{\tau} and all other slots of \code{\link{nl.fitt}} object which represent the variance function fitt information.}
    
    \item{\code{autcorr}:}{Object of class \code{"listorNULL"} of autocorrelated error.}
    
    \item{\code{autpar}:}{Object of class \code{"listorNULL"} of aprameters for autocorrelated error.}
    
    \item{\code{gresponse}:}{Object of class \code{"vectororMatrix"} generalized response, transformed response equal \eqn{R \times y}, for cholesky decomposition \eqn{R} of covariance matrix of error.}
    
    \item{\code{gpredictor}:}{Object of class \code{"vectororMatrix"} generalized predictor, transformed of predictor equal \eqn{R \times \eta(\theta)}, for cholesky decomposition \eqn{R} of covariance matrix of error.}
    
    Robust estimate slots:
    
    \item{\code{htheta}:}{Object of class \code{"vectororNULL"} optimized objective loss function is equal sum of rho function, with gradient and hessian as attribute.}
    
    \item{\code{rho}:}{Object of class \code{"vectororNULL"} computed robust \eqn{\rho} function, including gradient and hessian as attribute.}
    
    \item{\code{ri}:}{Object of class \code{"vectororNULL"} residuals equal predictor values minus predicted values, with gradient and hessian as attribute.}
    
    \item{\code{curvrob}:}{Object of class \code{"listorNULL"} robust Object of class \code{"listorNULL"} of PE and IE curvatures. Is not operational at the moment.}
    
    \item{\code{robform}:}{Object of class \code{"nl.formorNULL"}, robust \eqn{\rho} function of object type \code{"nl.form"}.}
    
    Nonlinear model estimates, inherited slots from \code{\link{nl.form}} object follows.
    
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
Class \code{"\link{nl.fitt.rob}"}, directly.
Class \code{"\link{nl.fitt}"}, by class "nl.fitt.rob", distance 2.
Class \code{"\link{nl.fitt.roborNULL}"}, by class "nl.fitt.rob", distance 2.
Class \code{"\link{nl.fittorNULL}"}, by class "nl.fitt.rob", distance 3.
}
\section{Methods}{
  \describe{
    \item{parInfer}{\code{signature(object = "nl.fitt")}: parameter inference function, calculate covariance matrix of parameters and their confidence interval. Usage: \code{parInfer(object,confidence = .95)}
    }
    \item{predictionI}{\code{signature(nlfited = "nl.fitt.gn")}: prediction interval. Usage: \code{predictionI(nlfited,confidence=.95,data=NULL)}, data is new data that will be predicting the values for them.}
    
    \item{residuals}{\code{signature(object = "nl.fitt.gn")}: residuals of fitt.}	 
  
    \item{atypicals}{\code{signature(nlfited = "nl.fitt")}: detect atypical points by calculating outlier detection measures. Usage: \code{atypicals(nlfited)}}
}
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  All information of a generalized nonlinear robust fited model are saved in \code{nl.fitt.rgn}, thus it can be large variable of informations. It is inheritance of \code{\link{nl.fitt.rob}}, and robust form of \code{\link{nl.fitt.gn}}. It include heterogeneous or autocorrelated fitt. The heteroscedastic fitt result stores in \code{hetro} slot, and autocorelation result stores in \code{correlation,autcor} slot. meanwhile the \code{vm, rm} include contains general of covariance and correlation matrix of both heteroscedastic and autocorrelated informations.
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{fittmethod}}, \code{\link{nl.fitt.rob}}, \code{\link{Fault}}, \code{\link{nl.form}}
}
\examples{
showClass("nl.fitt.rgn")
}
\keyword{classes}
\keyword{classes}
\keyword{robust fitt}
\keyword{least square}
\keyword{heteroscedastic variance}
\keyword{MM-estimates}
\keyword{Generalized least square}
\keyword{Generalized robust estimate}


