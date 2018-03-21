\name{db.Fault}
\alias{db.Fault}
\docType{data}
\title{Fault database
}
\description{Codes of errors used in "nlr" functions.
}
\usage{db.Fault}
\format{
  The format is:data.frame
 chr "db.Fault"
 \itemize{
  \item FL: (Fault Logic) is true if mesage is error and program terminate, False if message is warning and result might not be accurate.
  \item FN: (Fault  Number) is a code for the error mesage.
  \item FT: (Fault Text) is error text.
  \item FF: (Fault File) is the first origin of fault createor, but during the computation will represent the function that raised error.
 }
}
\details{db.Fault is the database of errors that might occure in all functions of the nlr package. Termination of the functions can be due to error or warning. If error happened no output is returned, if warning happened output returned but might not be reliable. Warning kight happen for example if number of iteration exceeded the maximum number of iteration. Then a sugestion will be displayed for remedy the problem. The returned variable is a "Fault" that display the mesage code, number, and file that error happened.
}
\source{Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al. 
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\examples{
data(db.Fault)   # load data
db.Fault        # variable name
}
\keyword{datasets}
