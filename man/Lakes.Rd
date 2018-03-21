\name{Lakes}
\alias{Lakes}
\docType{data}
\title{Lakes Data
}
\description{Lakes data represent a multivariate predictor nonlinear model.
}
\usage{data(Lakes)}
\format{
  The format is:data.frame
 chr "Lakes"
 \itemize{
    \item tn: mean annual total nitrogen concentration, response variable.
    \item nin: average influence nitrogen concentration, predictor variable one.
    \item tw: water retention time, predictor variable two.
 }
}
\details{Lakes data is collected from 29 lakes in Florida by United States Environmental Protection Agency (1978). Stromberg (1993) has identified observations 10 and 23 as outliers. The data presents the relationship between the mean annual total nitrogen concentration, TN, as the response variable and the average influence nitrogen concentration, NIN, and water retention time, TW, as predictors.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{United Stated Environmental Protection Agency (1978), "A Conpendium of Lake and Reservoir data Collected by the National Eutrophication Survey in Eastern, North Centeral and Southeastern United States," Working Paper #475, Corvallis Environmental Research Labratory, Corvallis, Oregon.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Stromberg, A. J. (1993). Computation of High Breakdown Nonlinear Regression Parameters, Journal of American Statistical Association 88(421): 237-244.
}
\examples{
data(Lakes)
Lakes
}
\keyword{datasets}
