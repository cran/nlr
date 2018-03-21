\name{nl.fitt.gn}
\Rdversion{1.1}
\docType{class}
\alias{nl.fitt.gn}
\alias{nl.fitt.gn-class}
\alias{atypicals,nl.fitt.gn-method}
\alias{parInfer,nl.fitt.gn-method}
\alias{predictionI,nl.fitt.gn-method}
\alias{recalc,nl.fitt.gn-method}
\alias{residuals,nl.fitt.gn-method}

\title{Class \code{"nl.fitt.gn"}}
\description{
  Object for Generalized fitt of nonlinear regression, inherited from \code{\link{nl.fitt}} object. Classic estimates of Heterogeneouity, Autocorrelated, or Weighted fitt will be saved in this object.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nl.fitt.gn",...)}. Include all slots of \code{nl.fitt} and another slots for saving heterosedastic variance fit or autocorrelated error parameters.
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
    So far was parent slot \code{\link{nl.fitt}}. Following slots are other that represent generalized nonlinear regression fitt.
    \item{\code{vm}:}{Object of class \code{"matrix"} of variance covariance matrix of error.}
    \item{\code{rm}:}{Object of class \code{"matrix"} of correlated error.}
    \item{\code{hetro}:}{Object of class \code{"nl.fittorNULL"}, include object \code{nl.fitt} of heteroscedastic error fit, or \code{NULL} for non hetroscedastic. It include parameter estimates of hetroscedastic variance \eqn{\tau} and all other slots of \code{\link{nl.fitt}} object which represent the variance function fitt information.}
    \item{\code{autcorr}:}{Object of class \code{"listorNULL"} of autocorrelated error.}
    \item{\code{autpar}:}{Object of class \code{"listorNULL"} of aprameters for autocorrelated error.}
    \item{\code{gresponse}:}{Object of class \code{"vectororMatrix"} generalized response, transformed response equal \eqn{R \times y}, for cholesky decomposition \eqn{R} of covariance matrix of error.}
    \item{\code{gpredictor}:}{Object of class \code{"vectororMatrix"} generalized predictor, transformed of predictor equal \eqn{R \times \eta(\theta)}, for cholesky decomposition \eqn{R} of covariance matrix of error.}
}
}
\section{Extends}{
  Class \code{"\link{nl.fitt}"}, directly.
  Class \code{"\link{nl.fittorNULL}"}, by class "nl.fitt", distance 2.
}
\section{Methods}{
  \describe{
    \item{atypicals}{\code{signature(nlfited = "nl.fitt.gn")}: Compute statistical measures to identify outliers.}
    \item{parInfer}{\code{signature(object = "nl.fitt.gn")}: parameter inference, covariance matrix of parameters.}
    \item{predictionI}{\code{signature(nlfited = "nl.fitt.gn")}: prediction interval.}
    \item{recalc}{\code{signature(object = "nl.fitt")}: generic function, recalculate the object with new arguments given in dots argument.}
    \item{residuals}{\code{signature(object = "nl.fitt.gn")}: residuals of fitt.}
	 }
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
  \code{nl.fitt.gn} inherit \code{nl.fitt} object. It include heterogeneous or autocorrelated fitt. The heteroscedastic fitt result stores in \code{hetro} slot, and autocorelation result stores in \code{correlation,autcor} slot. meanwhile the \code{vm, rm} include contains general of covariance and correlation matrix of both heteroscedastic and autocorrelated informations. See Seber and Wild (2003) for details in generalized Nonlinear model.
}
\seealso{
  \code{\link{nl.fitt}}, \code{\link{fittmethod}}.
}
\examples{
showClass("nl.fitt.gn")
}
\keyword{classes}
