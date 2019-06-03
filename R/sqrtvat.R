#*********************************************************************************
#**   function: sqrtvat, square root vairiance attribute.                       **
#**             transform variance to standard deviation                        **
#**             with all its gradient and hessian                               **
#**      variables:  vc: is (n*1) vector of some variance, transform to sqrt(vc)**
#**                  attr(vc,"gradient"), (n*P) gradient.                       **
#**                  attr(vc,"hessian")    (n,p,p) hessian                      **
#**      gradiet -> (1/2 \/vc) * grad(vc)                                       **
#**      hessian -> (1/2 \/vc ) * hesian(vc) - (1/4*sigma^3) grad(vc)T "%m3d"   **
#**                                                          grad(vc)           **
#*********************************************************************************
sqrtvat <- function(varcomp)
{
	vc <- as.numeric(varcomp$predictor)
	vcg <- attr(varcomp$predictor,"gradient")       # (n*p) gardient                     #
	vch <- attr(varcomp$predictor,"hessian")        # (n*p*p) hessian                    #
	stdv <- sqrt(vc)                                #  sigma i (standard deviation)      #
	stdvg <- (1.0/2*stdv) * vcg                     #  grad (sigmai)   gradient stdev    #
	
	vcgt <- t(vcg)
	gtg <- vcgt %m3d% vcg                           # vcgT * vc (n*p*p)                  #
	.temp2 <- 4.0 * vc * stdv                       # 4* sigma ^3                        #
	.expr2 <- (1.0 / .temp2) * gtg                  # 1/4* sigma^3 * gr T * grad         #
	.expr1 <- (1.0 / 2.0 * stdv) * vch              # 1/2*sigma  * hess                  #
	stdvh <- .expr1 - .expr2
	attr(stdv,"gradient") <- stdvg
	attr(stdv,"hessian") <- stdvh	
	return(stdv)
}
