\name{db.method}
\alias{db.method}
\docType{data}
\title{
  methods database
}
\description{
  Database for methods used in "nlr" package functions.
}
\usage{db.method}
\format{
  The format is:data.frame
 chr "db.method"
 \itemize{
  \item{\code{methodID}:}{ Object of class \code{"numeric"} code for the method.}
  
   \item{\code{method}:}{Object of class \code{"character"} name of the method used in estimation and other computation procedures.}
   
   \item{\code{detail}:}{Object of class \code{"character"} detail text description of the method.}
   
  \item{\code{methodBR}:}{Object of class \code{"numeric"} (method Branche) branch for the method used in iteration.}
  
  \item{\code{detailBR}:}{Object of class \code{"character"} detail of the branch.}
  
  \item{\code{subroutine}:}{Object of class \code{"character"} the function, subroutine, that the result constructed from.}
  
 }
}
\details{Each function use a method for estimation, and each method depending on numerical computation have a branch. Any output object have a "fittmethod" object which the values come from "db.method" database.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\examples{
data(db.method)        # load data
db.method             # variable name
}
\keyword{datasets}
