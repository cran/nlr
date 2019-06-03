\name{plotlist}
\alias{plotlist}
\title{
  plot a list of objects.}
\description{
  Plot each item of objects stored in a list. This is not extention of list object, just to help plot many objects stored in a list. Typically \code{\link{nlr}} fited objects such as \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}} has plot method and the result of many model fitts can be stored in a list, then \code{plotlist} function plot all of them.
}
\usage{
  plotlist(listobj,...)
}

\arguments{
  \item{listobj}{
    An \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}} obejct after a nonlinear model fitted before.
  }
  \item{...}{
    An \code{...}, argument passe to \code{plot} function.
  }
}
\details{
  plot items in a least one by one.
}
\value{
 Plot of objects.
    
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
  Used mostly for internal purposes.
}
\seealso{
 \code{\link{nl.fitt}}, \code{\link{nl.fitt.gn}}, \code{\link{nl.fitt.rob}}, \code{\link{nl.fitt.rgn}}
}
\examples{
nn<-length(nlrobj1)
d<-list(xr=nlr::Weights$Date, yr=nlr::Weights$Weight)
wmodel.nlrob <-list( NULL)
j<-1
for(i in 1:3){

  start<-as.list(nlrobj1[[i]]$par)
  #wmodel.nlrob[i] <-list( nlmest(nlrobj1[[i]],data=d,start=start,robfunc=nl.robfuncs[[1]]))
  a1 <- nlr(nlrobj1[[i]],data=d,robustform = "hampel")
  if (! is.Fault(a1)){
    wmodel.nlrob[j]<-a1
    j<-j+1
  }
}
plotlist(wmodel.nlrob)
}
\keyword{plot}
\keyword{plotlist}
