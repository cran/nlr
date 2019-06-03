#singularCase=2 means compute   zeta <- hsh + lambda * diag(abs(diag(hsh)))
# in singular case, like optim.NLM

nlr.control <-
  function(maxiter = 50,
           tolerance = 0.0001,
           minscale = 0.001,
           trace = F,
           minlanda = 1e-16,
           derivfree = F,
           robscale = T,
           algorithm = c("Levenberg-Marquardt", "Nelder-Mead", "Gauss Newton"),
           method = c("default",
                      "RME",
                      "CME",
                      "CLSME",
                      "RGME",
                      "WME",
                      "MLE",
                      "OLS",
                      "TS",
                      "RTS",
                      "lms"),
           initials = c("manuall", "lms", "OLS", "quantile"),
           history = F,
           length.out = NULL,
           singlePlot = F,
           singularCase = 1,
           JacobianLeverage = c("default", "classic", "robust")
           ) {
    algorithm <- match.arg(algorithm)
    method <- match.arg(method)
    initials <- match.arg(initials)
    JacobianLeverage<-match.arg(JacobianLeverage)
    list(
      maxiter = maxiter,
      tolerance = tolerance,
      minscale = minscale,
      minlanda = minlanda,
      trace = trace,
      derivfree = derivfree,
      robscale = robscale,
      method = method,
      algorithm = algorithm,
      initials = initials,
      history = history,
      length.out = length.out,
      singlePlot = singlePlot,
      singularCase = singularCase,
      JacobianLeverage=JacobianLeverage
    )
  }
