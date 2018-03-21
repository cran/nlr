\name{convfkt2nlform}
\alias{convfkt2nlform}

\title{
  Convert \code{fktlist} objects to \code{nl.form}.
}
\description{
Convert \code{fktlist} objects defined by Bunke et.al to \code{nl.form}. It calculate gradient and hessian using derive3 function.
}
\usage{
convfkt2nlform(fktlistex, namesdata = NULL)
}
\arguments{
  \item{fktlistex}{ fktlistex object include fkt, par, p, ccode, lambda, case, inv, name, defined by Bunke et al (1998)
}
  \item{namesdata}{
  Name of data, is not functional now, used for feature development.
}
}
\details{
\code{fktlist} defined by Bunke et al (1998) \code{convfkt2nlform} function convert this object to \code{nl.form} object.
}
\value{
\code{nl.fomr} object with gradient and hessian attributes.
}
\references{
 Bunke, O., Droge, B., Polzehl, J. Splus tools for model selection in nonlinear regression (1998) Computational Statistics, 13 (2), pp. 257-281.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
Derivatives must not be included in \code{fktlistext} object will be added as attributes to response or independent variables.
}

\seealso{
See Also as \code{\link{convexpr2nlform}}, \code{\link{nl.form}}.
}
\examples{
## fktlistex created by Bunke et.al
convfkt2nlform(fktlistex1[[1]])
}
\keyword{ nl.form }
\keyword{ Gradient }
\keyword{ Hessian }
