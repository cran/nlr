\name{indifinv}
\alias{indifinv}
\title{Indefinite Inverse of matrix.
}
\description{Compute the inverse of matrix using spectoral decomposition, using eigenvalues and eigen vectors of matrix.
}
\usage{
indifinv(mtrx, stp = T, symmetric = all(mtrx == t(mtrx)))
}
\arguments{
  \item{mtrx}{square matrix to compute the inverse.
}
  \item{stp}{if stp=T when error happened stop running program, if stp=F, does not stop program but return back Fault object.
}
  \item{symmetric}{symmetric}{Used for computing eigenvalues, if symmetric=T the matrix is symetric, if symmetric=F the matrix is not symetric.
}
}
\details{\code{eiginv} function compute the inverse of matrix using spectoral decomposition 
  \deqn{
   A_{k \times k}=\textbf{P} \Lambda \textbf{P}'
  }
  where
  \deqn{
  	\textbf{P}=[e_1,\dots,e_k]
  }
  \deqn{
  	\Lambda=diag(\lambda_i)
  }
  in which \eqn{\lambda_i} is eigenvalues of matrix A coresponding to eigenvector       \eqn{e_i}. Then the inverse is:
  \deqn{
	  A^{-1}=\textbf{P} \Lambda^{-1} \textbf{P}'
  }
}
\value{
  If matrix eigenvalues are not zero return the inverse of matrix, otherwise returns \code{\link{Fault}} object with fault number=9, means the matrix is not positive definit.
}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Jan 2010. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{
  \code{\link{eiginv}} return back error if negative eigen values occured, means singular matrix. But \code{indifinv} raise error if infinite or null values occured. The eigenvalues decomposition is used as \code{\link{eiginv}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{eiginv}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 a1=matrix(c(1,2,3,4,5,6,7,8,9),nrow=3)
 indifinv(a1)
}
\keyword{Matrix Inverse}
\keyword{Indifinit Inverse}
