\name{parameter.names}
\alias{parameter.names}
\title{
Support for Functions \code{nlr()}
}
\description{
This is support for the functions ms() and nls(). It is not intended to be called directly by users. 
}
\usage{
parameter.names(formula, data)
}
\arguments{
  \item{formula}{
    Expression include variables.
  }
  \item{data}{
    data might have parameters attribute.
}
}
\details{
  In several functions, If start is not given the name of data and variables have to be gathered from formula and environment variables, which use this function appropriately.
}
\value{
  list of parameters.
}
\references{
Venables, W.N., and Ripley, B.D. (1999). Modern Applied Statistics with S-PLUS. New York: Springer-Verlag
}
\author{
Venables, W.N., and Ripley
}
\note{
For internal use, muight not be called by user.
}
\seealso{
\code{\link{nlmest.NLM}}
}
\examples{
## The function is currently defined as
"parameter.names"
}
\keyword{ parameters }
\keyword{ expression}% __ONLY ONE__ keyword per line
