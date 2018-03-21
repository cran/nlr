\name{mplot}
\alias{mplot}
\title{
  Multiple Plot
}
\description{
  plot multiple models stored in \code{nl.fitt} and its child objects, in same graph.
}
\usage{
  mplot(mlist, case = 1, length.out = NULL, ...)
}

\arguments{
  \item{mlist}{
    list of object  models
  }
  \item{case}{
    case=1, common x, case=2, different x
  }
  \item{length.out}{
    length of predictor to be increamented, if not given the original predictor data will be used.
  }
  \item{\dots}{
    extra option submit to plot.
}
}
\details{
  If different methods are used to estimate a model this function can be used to plot them all in same graph. Multiple output of estimates in \code{nlr} package should be stored in a list, then \code{mplot} plot data and predicted values for all methods over the data.
}
\value{
Plot graph.
}
\references{
Rizo ML 2008 Statistical Computing with R The R Series. Chapman & Hall/CRC The R Series.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{www.riazoshams.com/nlr}
}
\note{
  A more common situation is when case=1 then common x-axis will be considered for all fited objects. That is the case when sommon data used with different method of fitt and once want to compare different methods.
}

\seealso{
  \code{nl.fitt}
}
\examples{
datalist<-list(xr=trade.ir[,1],yr=trade.ir[,2])
a1<- nlr(nlrobj5[[18]],data=datalist,
control=nlr.control(trace=TRUE,derivfree = FALSE,tolerance=1e-4,singularCase=1,maxiter = 8000))    
a2<- nlr(nlrobj5[[18]],data=datalist,
control=nlr.control(trace=TRUE,method="OLS",
derivfree = FALSE,tolerance=1e-4,maxiter = 4000))    # nelder mead, MM, selfstart

mlist=list(a1,a2)
mplot(mlist)
}
\keyword{Multiple plot}

