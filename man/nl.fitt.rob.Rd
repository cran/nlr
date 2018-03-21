\name{nl.fitt.rob-class}
\Rdversion{1.1}
\docType{class}
\alias{nl.fitt.rob}
\alias{nl.fitt.rob-class}
\alias{dlev}
\alias{dlev,nl.fitt.rob-method}
\alias{JacobianLeverage,nl.fitt.rob-method}
\alias{parInfer,nl.fitt.rob-method}
\alias{plot,nl.fitt.rob-method}
\alias{plot,nl.fitt.rob,missing-method}
\alias{plotlms}
\alias{predictionI,nl.fitt.rob-method}
\alias{recalc,nl.fitt.rob-method}

\title{Class \code{"nl.fitt.rob"}}
\description{
  Object of robust estimates of nonlinear regression model.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nl.fitt.rob", ...)}.
}
\section{Slots}{
  \describe{
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
Class \code{"\link{nl.fitt}"}, directly.
Class \code{"\link{nl.fitt.roborNULL}"}, directly.
Class \code{"\link{nl.fittorNULL}"}, by class "nl.fitt", distance 2.
}
\section{Methods}{
  \describe{
        \item{dlev}{\code{signature(nlfited = "nl.fitt.rob")}: DLEV Difference in LEverage measure.}
        \item{JacobianLeverage}{\code{signature(nlfited = "nl.fitt.rob")}: Jacobian-Leverage for nonlinear regression. Usage \code{JacobianLeverage(nlfited = "nl.fitt.rob")}}
        
    \item{parInfer}{\code{signature(object = "nl.fitt")}: parameter inference function, calculate covariance matrix of parameters and their confidence interval. Usage: \code{parInfer(object,confidence = .95)}
    }
    
    \item{plot}{\code{signature(x = "nl.fitt", y = "missing",control=nlr.control(history=F,length.out=NULL,singlePlot=F),...)}: generic function extended to \code{nl.fitt} object. Plot the object. Usage. \code{plot(x,y="missing",control=nlr.control()}. If history is TRUE the convergence of fitt will be ploted.length.out is length of incremented p[redictor to acheive smooter curve. singlePlot=F plot the model and residuals in two collumn.
    If the estimate be Least MEdian Square, the \code{plotlms} function is used to plot the object.
    }

    \item{predictionI}{\code{signature(nlfited = "nl.fitt.gn")}: prediction interval. Usage: \code{predictionI(nlfited,confidence=.95,data=NULL)}, data is new data that will be predicting the values for them.}
    \item{recalc}{\code{signature(object = "nl.fitt.rob")}: recalculate the original call of the fited model by some extra options. It is created for usage in \code{atyp} function. }
  }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons..
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  All information of a nonlinear robust fited model are saved in \code{nl.fitt.rob}, thus it can be large variable of informations. The generalized form \code{\link{nl.fitt.rgn}} of a fitt is children of this object. Typically it used to save robust MM-estimation method. But it is extensively used to save other fited objects such as hetroscedastic variance parameter fitts.
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{fittmethod}}, \code{\link{nl.fitt.rgn}}, \code{\link{Fault}}, \code{\link{nl.form}}
}
\examples{
showClass("nl.fitt.rob")
}
\keyword{classes}
\keyword{robust fitt}
\keyword{least square}
\keyword{heteroscedastic variance}
\keyword{MM-estimates}

