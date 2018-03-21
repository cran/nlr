\name{db.methodBR}
\alias{db.methodBR}
\docType{data}
\title{method branches database
%%   ~~ data name/kind ... ~~
}
\description{Sub method used in procedures in functions.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{db.methodBR}
\format{
  The format is: data.frame
 chr "db.methodBR"
 \itemize{
    \item methodBR: method branch.
    \item detailBR: detail and long description of the method.
 }
}
\details{Any method used in functions might have some sub branch depending the procedure it used. The codes and detail is from "db.methodBR" database. The result objects from function have "fittmethod" sub object that save methods and sub methods used in the function.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{Robust Nonlinear Regression, Theories and Methods with Practical Guides for R Packages. Riazoshams et al.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\examples{
data(db.methodBR)       # load data
db.methodBR            # variable name
## maybe str(db.methodBR) ; plot(db.methodBR) ...
}
\keyword{datasets}
