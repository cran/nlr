\name{carbon}
\alias{carbon}
\docType{data}
\title{Carbon Dioxide data.
}
\description{Carbon dioxide traped in iceberg during history.
}
\usage{data(carbon)}
\format{
  The format is: data.frame
 chr "carbon"
 \itemize{
    \item year: year of gas traped in iceberg.
    \item co2: measured Carbon Dioxide.
 }
}
\details{UNEP (1989) presented the Methane Gas and Carbon Dioxide Gas collected from the Gas trapped in icebergs in south  pole from 8000 years ago. 
}
\source{
  UNEP (1989), Environmental data report / prepared for UNEP by the GEMS Monitoring and Assessment Research Centre, London, UK, in co-operation with the World Resources Institute, Washington, D.C.
}
\references{
  Riazoshams, H., Miri, H., (2013) Application of Robust Nonlinear Regression, case study for modeling the greenhouse gases, Methane and Carbon Dioxide concentration in atmosphere. International Conference on Mathematical Science and Statistics (ICMSS 2013), Kula Lumpur, Malaysia.
}
\examples{
data(carbon)
 carbon$year
 carbon$co2
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
plot(carbon.ols,control=nlr.control(history=TRUE))
}
\keyword{datasets}
