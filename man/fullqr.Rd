\name{fullqr}
\alias{fullqr}
\title{
full rank QR decomposition.}
\description{
Is actual QR decomposition matrix, created for compatibility to \code{nlr} package functions.
}
\usage{
fullqr(x)
}

\arguments{
  \item{x}{
A Square matrix.
}
}
\details{
Compute QR matrices, and provide explicit matrix form of lower triangular and upper triangular matrices.
}
\value{
  List of standard matrix form to be used in internal functions in \code{nlr}.
  \itemize{
  \item{q}{\code{q} matrix.}
  \item{r}{\code{r} matrix.}
  \item{q2}{\code{q2} is from partition of Q=[q1|q2] matrix.}
  \item{r1} {\code{r1} is from partition of R=[R1/0] matrix.}
  \item{rinv}{is generalized inverse of \code{r1}.}
  \item{ar}{QR decomposition result from \code{qr} function.}
}
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Jan 2010.
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
Used mostly for internal purposes.
}
\seealso{
  \code{\link{nlsqr}}
}
\examples{
 a1=matrix(c(1,2,3,4,5,6,7,8,9),nrow=3)
 fullqr(a1)
}
\keyword{QR Decomposition}
\keyword{Full Rank QR decomposition}
