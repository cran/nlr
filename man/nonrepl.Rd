\name{nonrepl}
\alias{nonrepl}
\title{
  Sample variance of response.
}
\description{
  Compute the sample standard deviation of response data$y at repeated predictor data$x.
}
\usage{
nonrepl(data)
}
\arguments{
  \item{data}{
    list of response data$y and predictor data$x.
}
}
\details{
  If predictor \eqn{x(i)} repeated \eqn{n_i} times, the sample variance of response is compute for \eqn{y(i,j), j=1,...,n_i}.
}
\value{
  list of result include standard deviations and information about repeated data as:
  \item{x}{ordered data predictor.}
  \item{y}{ordered data response by predictor.}
  \item{xk}{nonrepeated data of x's.}
  \item{ni}{number of repeated of each xk's.}
  \item{xm}{position of each xk's in ordered x.}
  \item{k}{length of non repeated data, xk.}
  \item{xo}{}
  \item{yq}{}
  \item{ys}{}
}
\references{
  This function provided by
  Bunke, O., Droge, B., Polzehl, J. Splus tools for model selection in nonlinear regression (1998) Computational Statistics, 13 (2), pp. 257-281.
}
\author{
  Bunke et al. 1998.
}

\seealso{
  \code{\link{zvalues}}, \code{\link{rzvalues}}
}
\examples{
## The function is currently defined as
"nonrepl"
}
\keyword{repeated measures}

