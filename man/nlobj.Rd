\name{nlrobj1}
\alias{nlrobj1}
\alias{nlrobj3}
\alias{nlrobj4}
\alias{nlrobj5}
\alias{nlrobj6}
\alias{nlrobj7}
\alias{fktlistex1}
\alias{fktlistex3}
\alias{fktlistex4}
\alias{fktlistex5}
\alias{fktlistex6}
\alias{fktlistex7}
\alias{fktmodels1}
\alias{fktmodels3}
\alias{fktmodels4}
\alias{fktmodels5}
\alias{fktmodels6}
\alias{fktmodels7}
\docType{data}
\title{
  Nonlinear model objects
}
\description{
  List of \code{\link{nl.form}} objects of nonlinear regression models.
  \itemize{
  \item{nlrobj1:}{ has 16 model.}
  \item{nlrobj3:}{ has 18 model.}
  \item{nlrobj4:}{ has 11 model.}
  \item{nlrobj5:}{ has 19 model.}
  \item{nlrobj6:}{ has 7 model.}
  \item{nlrobj7:}{ has 23 model.}
  }
}
\usage{
  nlrobj1
  }
\format{
  The format is: chr "nlrobj1"
  
  Elements are list of model, each cell is \code{\link{nl.form}} object.
}
\details{
by index can access any of the nodels, for example nlrobj1[[14]] with name nlrobj1[[14]]$name="Logistic without intercept" is logestic model to fitt the chicken \code{\link{Weights}} data.
}
\source{
  ("nlrobj1", "nlrobj3", "nlrobj4", "nlrobj5") are generalization of objects developed by Bunke et al. 1998.
}

\references{
Bunke, O., Droge, B., Polzehl, J. Splus tools for model selection in nonlinear regression (1998) Computational Statistics, 13 (2), pp. 257-281.
}
\note{
 The variables with prefix \code{fkt} are created by Bunke et al. They will not be used in \code{nlr} directly or by user.
}
\seealso{
  \code{\link{nl.form}}, \code{\link{Weights}}
}
\examples{
data(nlrobj1)
nlrobj1
}
\keyword{datasets}
