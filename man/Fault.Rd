\name{Fault}
\Rdversion{1.1}
\docType{class}
\alias{Fault}
\alias{Fault-class}
\alias{$,Fault-method}

\title{Class \code{"Fault"}}
\description{
  Error object returns from all functions in \code{nlr}. It may return no fault or warning or error message.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Fault",FL = F,FN=0,FT=NULL,FF=NULL,pnt=F)}.
  
  if \code{pnt=T} then it display the created object result.
}
\section{Slots}{
  \describe{
    \item{\code{FL}:}{Object of class \code{"logical"} Fault Logic if is TRUE error occured otherwise no error occured.}
    \item{\code{FN}:}{Object of class \code{"numeric"} Fault Number, code for error, see notes bellow. All codes can be seen in \code{\link{db.Fault}}
    }
    \item{\code{FT}:}{Object of class \code{"character"} Fault Text, explain the error.}
    \item{\code{FF}:}{Object of class \code{"character"} Fault File, the function that error raised in.}
  }
}
\section{methods}{
  \describe{
    \item{$}{\code{signature(x = "Fault")}: return slots.}
	 }
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  \code{\link{db.Fault}} variable include error codes and expression used in \code{Fault}.
}
\seealso{
  \code{\link{db.Fault}}
}
\examples{
showClass("Fault")
db.Fault        # Error codes
Fault(FN=1)
}

