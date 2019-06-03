#*****************************************************************************************************
individ <- function(hessian)
{
  n <- length(hessian[, 1, 1])
  p <- length(hessian[1, 1, ])
  sp <- matrix(hessian[, 1, 1], ncol=1)
  for (i in 2:p)
  {
    sp <- cbind(sp, hessian[, i, 1:i])
  }
  return(sp)
}
