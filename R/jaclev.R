
#+------------------------------------------------------------------------------+
#|    jaclev: Compute the Jacobian Leverage, generalized for nonlinear case.    |
#|        by Roy, Cook, 1992, JASA.                                             |   
#|    Inputs:                                                                   |
#|        gradient, (n*p) gradient.                                             |
#|        hessian, (n*p*p) hessian.                                             |
#|        rsd, (n*1) residual vector.                                           |
#+------------------------------------------------------------------------------+

jaclev <- function(gradient,hessian,rsd){
#for(i in 1:dim(hessian)[2])
	prva <-   prodVA(hessian,rsd)                ##  [e] [w]   p*p
	.expr1 <- t(gradient) %*% gradient           ##  vT * v    p*p
	.expr2 <- .expr1 - prva                      ##  vT v - [e] [w]    p*p
	.expr3 <- indifinv(.expr2,stp=F,symmetric=T)   ##  (vT v - [e] [w])^-1   p*p
	if(class(.expr3)=="Fault") return(.expr3)
	.temp1 <- gradient %*% .expr3                ##  v * (vT v - [e] [w])^-1    n*p
	JacLeverage <- .temp1 %*% t(gradient)        ##  v * (vT v - [e] [w])^-1 vT  n*n
	return(JacLeverage)
}