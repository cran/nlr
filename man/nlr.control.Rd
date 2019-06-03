\name{nlr.control}
\alias{nlr.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{list of nlr package controls.
}
\description{create nlr.control options for using in diferent part of the package.
}
\usage{
nlr.control(maxiter = 50, tolerance = 0.0001, minscale = 0.001, trace = F,
minlanda=1e-16,derivfree=F,robscale=T,
algorithm = c("Levenberg-Marquardt", "Nelder-Mead","Gauss Newton"),
method=c("default","RME","CME","CLSME","RGME","WME","MLE","OLS","TS","RTS","lms"),
initials=c("manuall","lms","OLS","quantile"),history=F,length.out=NULL,singlePlot=F,
singularCase=1,
JacobianLeverage = c("default", "classic", "robust"))
}
\arguments{
  \item{maxiter}{maximum number of iteration. To be used optimization or fitt procedures.
}
  \item{tolerance}{Tolerance of convergence.
}
  \item{minscale}{Minimum of scale value.
}
  \item{trace}{Default is False, if True draw the graph of convergence in iterations, depends on availability for that cases.
}
  \item{minlanda}{Minimum of landa coeficient values in stepest descend or levenberg marquardt. When the step length does not reduce the objective function this value decrease. 
}
  \item{derivfree}{Default value is False. If True the implicitly program use derivative free optimization methods based on Nelder-Mead method. This is very slow convergence in that case the number of maximum iteration should be increased, at least maxiter=500 is sugested. nlr is derivative based but in the cases where problem happenes in derivative computation Nelder-Mead will be used.
}
  \item{robscale}{Default is True, when False the clasic standard error of residuals will be used in computations. It is not recomneded to use this option because reduce the robustness of estimators.
}
  \item{algorithm}{of optimization method, default is "Levenberg-Marquardt", is robust in computation when outlier happense. "Nelder-Mead" is derivative free method, and "Gauss Newton" method is used for fast computation but might face with some singularity in hessian when outlier happense.
  }
  \item{singularCase}{Select how to solve the singular gradient matrix case in MM-estimate procedure of levenberg marquard method. The efault value equal 1 add up a value to diagonal elements, and 2 add up proportion to size of diagonal values. They are same but Generally case 1 works in most of cases, but in case of divergence or singularity problem can use 2.}
  \item{JacobianLeverage}{character name of jacobian leverage to be used incomputation.
    \itemize{
      \item{"default":}{default value assigned by any function}
      \item{"classic":}{cllasic nonrobust value}
      \item{"robust":}{robust jacobian leverage value}
    }
  }
  \item{method}{of computation using several type of estimators.
    \itemize{
      \item {"RME":} {Heteroscedastic error Robust Multi Stage Estimate.}
      \item {"CME":} {Heteroscedastic error Classic Multi Stage Estimate.}
      \item {"CLSME":} {Heteroscedastic error Classic Least Square based Estimate.}
      \item {"RGME":} {Heteroscedastic error Robust Generalized Estimate.}
      \item {"WME":} {Heteroscedastic error Weighted M-Estimate, See Lim 2010.}
      \item {"MLE":} {Maximum Likelihood Estimate.}
      \item {"OLS":} {Ordinary Least square, constant and uncorrelated error.}
      \item {"TS":} {Autocorelated error, classic Two Stage Estimate.}
      \item {"RTS":} {Autocorelated error, Robust Two Stage Estimate.}
      \item {"lms":} {Least median of squared residuals estimate. (Non efficient)}
      }
  }
  \item{initials}{ is used to define initializing parameters using a specific estimatro. Might be used when initial values are dificult to find manually, or getInitial function.
    \itemize{
      \item {"manuall"} {default value is manually provided initial values by user.}
      \item {"lms"} {robust high breakdown point least median of squares of errors estimator.}
      \item {"ols"} {ordinary least square estimator.}
      \item {"quantile"} {robust least quantile regression.}
    }
    }
  \item{length.out}{Length for incrementing independent variable to be used in ploting commands to have a more smoother curve.}
  \item{history}{Default=F, Used in plot command, to draw the history of convergence, default is F, the TRUE value force the function to draw the history.}
  \item{singlePlot}{Default=F, By default the plot function draw the fited model (in one dimentional case), and residuals in second collumn, singlePlot=T causes the plot to draw the two graph in two windows.
  }
}
\details{nlr.control function create nlr.control variable and mostly used as control argument to the functions. It define convergence and many other parameters depends on the function operation.
%%  ~~ If necessary, more details than the description above ~~
}
\value{list of controls.

}
\references{
   Riazoshams H, Midi H, and Ghilagaber G, 2018,. Robust Nonlinear Regression, with 
Application using R, Joh Wiley and Sons.
}
\author{
  Hossein Riazoshams, Apr 2014. 
  Email: \email{riazihosein@gmail.com} 
  URL \url{http://www.riazoshams.com/nlr/}
}

\note{"nlr" package is derivative based, the default value for algorithm for optimization is derivative based methods such as Levenberg-Marquardt, then nlr check use optional user request.
Beside some options for ploting and further more can be defined by \code{nlr.control}.
This control passe to functions for controling computations, by users.
}
\seealso{
    \code{\link{nlr}}
  }
\examples{
## The function is currently defined as
"nlr.control"
nlr.control()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Tolerance}
\keyword{OLS-Estimate}% __ONLY ONE__ keyword per line
\keyword{MM-Estiamte}
\keyword{Levenberg-Marquardt}
\keyword{Nelder-Mead}
\keyword{Gauss Newton}