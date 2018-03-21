\name{Tumor}
\alias{Tumor}
\docType{data}
\title{
Tumor metastasis data.
}
\description{
The response variable is the fraction of breast cancer patients with metastases and the predictor variable is the tumor size. 
}
\usage{data("Tumor")}
\format{
  The format is:
 chr "Tumor"
   \itemize{
    \item {tumorSize}  {predictor variable is the tumor size.}
    \item {metastasized} {response variable is the fraction of breast cancer patients with metastases.}
}
}
\details{
The data consist of 12 observations. The response variable is the fraction of breast cancer patients with metastases and the predictor variable is the tumor size. 
}
\source{
Michaelson JS, Halpern E, Kopans D. Breast cancer: Computer simulation methodfor estimating optimal intervals for screening. Radiology. 1999;21:551-560.
}
\references{
Michaelson JS, Halpern E, Kopans D. Breast cancer: Computer simulation methodfor estimating optimal intervals for screening. Radiology. 1999;215:51-560.
}
\examples{
data(Tumor)
Tumor
}
\keyword{datasets}
