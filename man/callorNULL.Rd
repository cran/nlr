\name{callorNULL-class}
\Rdversion{1.1}
\docType{class}
\alias{callorNULL-class}
\alias{characterorNULL-class}
\alias{expressionorNULL-class}
\alias{fittmethodorNULL-class}
\alias{functionorNULL-class}
\alias{integerorNULL-class}
\alias{listorNULL-class}
\alias{logicalorNULL-class}
\alias{matrixororNULL-class}
\alias{nl.fitt.roborNULL-class}
\alias{nl.fittorNULL-class}
\alias{nl.formorNULL-class}
\alias{numericorNULL-class}
\alias{vectororMatrix-class}
\alias{vectororNULL-class}
\title{Class \code{"or classes"}}
\description{
A set of or classes for compatibility purposes. It is used to be able to identify null values.
\itemize{
  \item{callorNULL}{call or NULL class.}
  \item{characterorNULL} {character or null.}
  \item{expressionorNULL} {expression or null.}
  \item{fittmethodorNULL} {fittmethod or null.}
  \item{functionorNULL} {function or null.}
  \item{integerorNULL} {integer or null.}
  \item{listorNULL} {list or null.}
  \item{logicalorNULL} {logical or null.}
  \item{matrixororNULL} {matrix or null.}
  \item{nl.fitt.roborNULL} {nl.fitt.rob or null.}
  \item{fittorNULL} {fitt or null.}
  \item{nl.formorNULL} {nl.form or null.}
  \item{nl.numericorNULL} {numeric or null.}
  \item{vectororMatrix} {vector or Matrix.}
  \item{vectororNULL} {vector or null.}
}
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{methods}{
No methods defined with class "callorNULL" in the signature.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, 2013. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
or classes are created in \code{nlr} for compatibility purpose with \code{splus}.
}
\seealso{
\code{\link{nlr}}
}
\examples{
setClassUnion("nl.fittorNULL", c("nl.fitt", "NULL"))
}
\keyword{classes}
