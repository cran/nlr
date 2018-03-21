\name{evald}
\alias{evald}
\title{\code{eval} function.
}
\description{Implemented form of \code{eval} function with dots argument.
}
\usage{
evald(expr, envir = parent.frame(), 
enclos = if (is.list(envir) || is.pairlist(envir)) parent.frame() else baseenv(), ...)
}
\arguments{
  \item{expr}{an object to be evaluated. 
}
  \item{envir}{  
the environment in which expr is to be evaluated. May also be NULL, a list, a data frame, a pairlist or an integer as specified to sys.call.
}
  \item{enclos}{Relevant when envir is a (pair)list or a data frame. Specifies the enclosure, i.e., where R looks for objects not found in envir. This can be NULL (interpreted as the base package environment, baseenv()) or an environment.
}
  \item{\dots}{other arguments pass to expr.
}
}
\details{this function works exactly similar to \code{eval} except accept \dots argument, it is a compatibility adjustment to SPLUS.
}
\value{The result of evaluating the object: for an expression vector this is the result of evaluating the last element.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Apr 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}

\note{
  This is a specialisd for \code{nlr} package to imitate \code{eval} function for special uses, is not intended to be used by user.
}
\seealso{
  \code{\link{eval}}, \code{\link{nl.form}}
}
\examples{
## The function is currently defined as
"evald"
}
\keyword{ expression }
\keyword{ environment }% __ONLY ONE__ keyword per line
