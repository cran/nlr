\name{convexpr2nlform}
\alias{convexpr2nlform}
\title{
  Convert \code{expression} to \code{nl.form} 
}
\description{
  Convert two sided (or one sided) expression formula to \code{nl.form} object using \code{derive3} from \code{MASS} library.
}
\usage{
convexpr2nlform(form, namesdata=NULL, start, inv = NULL, name="User Defined",...)
}
\arguments{
  \item{form}{
    Must be one sided expression (defined by ~formula) or two sided (response~predictor), nonlinear regression function, include parameters, response and predictor variables.
}
  \item{namesdata}{
    optional character vector of name of data include independent and possibly dependent in two sided fomula. 
}
  \item{start}{
    list of parameters, for which the gradinet and hessian will be computed.
}
  \item{name}{A character name for the model}
  \item{inv}{inverse of the nonlinear functin model}
  \item{...}{Ane extra argument pass to \code{nl.form}}
}
\details{
  \code{nlr} package is gradient based algorithm, is based \code{\link{nl.form}} object in which gradient and hessian is available. If a nonlinear regression model formula is one sided or two sided formula and its gradient and hessian exist, the \code{convexpr2nlform} convert it to \code{\link{nl.form}} object by calling \code{derive3} from \code{MASS} library. Although the existence of derivative is strong assumption but using advance programs can acheive high precision computing.
}
\value{
  \code{\link{nl.form}} object of the nonlinear regression function.
  \item{formula:}{ formula one sided or two sided with gradinet and hessian as attribute.}
  \item{formtype:}{="formula"}
  \item{p:}{=length(start) is number of parameters.}
  \item{name:}{="User Defined"}
  \item{par:}{=start parameters.}
  \item{dependent}{ character vector of name of dependent variable.}
  \item{independent:}{character vector of name of independent variable.}
  \item{origin:}{=form}
}
\references{
Rizo ML 2008 Statistical Computing with R The R Series. Chapman & Hall/CRC The R Series.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  If the derivatives does not exist in  \code{\link{nlr}} function eplicitly the derivative option must set to derivative free.
  The \code{namesdata} is not functional in this version, implemented for further development. The name of parameters will be constructed from \code{start} arguments and the name of independent and dependent variables will be derived from the rest of variables embeded in the \code{form} expression.
}

\seealso{
\code{\link{nl.form}}, \code{\link{nlr}} 
}
\examples{
## The function is currently defined as
nlf=convexpr2nlform(yr ~ (a)*(exp(-b*xr)-exp(-c*xr)), start = list(a=.05,b=4.39,c=21.6))
nlf
}
\keyword{gradient}
\keyword{hessian}
\keyword{expression}
\keyword{nonlinear regression model}


