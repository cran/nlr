\name{curvature}
\alias{curvature}
\title{Calculate (IE) Intrinsic curvature and (PE) Parameter curvature.
}
\description{IE and PE are measures to identify the linear approximation of nonlinear model is appropriate or no. This function may not be called explicitly by user.
}
\usage{
curvature(gradient, hessian, sigma)
}
\arguments{
  \item{gradient}{n by p gradient of fited model.
}
  \item{hessian}{n by p by p array of hessian for the nonlinear model.
}
  \item{sigma}{estimated standard deviation.
}
}
\details{Gauss Newton method of estimation is based on linear approximation to nonlinear model. The linear approximation to function might not be appropriate. PE and IE is used to identify the parameter effect and intrinsic effect of model. Big values represent the linear approximation to nonlinear model is not correct.
}
\value{
  List of curvature values. 
  \itemize{
    \item pe: Parameter Effect curvature.
    \item int: Intrinsic effect curvature.
    \item a: A matrix.
    \item cutf: cut of point \deqn{1/\sqrt{F(.95,p,n-p)}}, if PE or IE be biger than cut of point then either of them has larg curvature.
  }
}
\references{Bates, D.M., and Watts, D. G. (1980). Relative curvature measures of nonlinearity, J. R. statistic. Ser. B 42: 1-25.
}
\author{
  Hossein Riazoshams, May 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}
\note{
  curvature is a model checking tool. From the OLS estimate output included carvature that can be accessed by \code{curvature} slot of the output object, therefore do not need to be called explicitly by user.
}
\examples{
 crbdt<-list(xr=nlr::carbon$year,yr=nlr::carbon$co2)
 ScalExp<- convexpr2nlform(yr ~ p1 + exp(-(p2 - p3 * xr)),
                        selfStart=function(data){
                            y1 <-as.double(data$yr)
                            p1<-min(y1)
                            y<-log(y1-p1+10*.Machine$double.eps)
                            x<-as.double(data$xr)
                            b1<-lm(y~x)
                            p2<- -b1$coefficients[1]
                            p3<- b1$coefficients[2]
                            return(list(p1=p1,p2=p2,p3=p3))
                        },
                        name="Scaled Exp convex",
                        start=list(p1=700,p2=21,p3=0.01)
)
 carbon.ols <- nlr(formula=ScalExp, data=crbdt, 
                  control=nlr.control(method="OLS"))
carbon.ols$curvature                  
}
\keyword{Curvature}
