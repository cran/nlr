#*************************************************************************************
#***                            Compute Curvatures, pe nd IN.                      ***
#***                                                                               ***  
#***                                                                               ***  
#***                                                                               ***  
#*************************************************************************************
curvature <- function(gradient,hessian, sigma)
{
	p <- ncol(gradient)
	n <- nrow(gradient)
	hs <- individ(hessian)
	dmat <- cbind(gradient,hs)
	qrd <- qr(dmat)
	if(qrd$rank < p) return(NULL)
	qd <- qr.Q(qrd)
	rd <- qr.R(qrd)
	rnk <- qrd$rank
	q1d <- qd[, 1:rnk]
	if(rnk==p) return(NULL)
	a <- array(0, c(rnk, p, p))
	for (j in 1:p)
	{
		a[, , j] <- crossprod(q1d, hessian[, ,j])
	}
	
	a <- aperm(a, c(2, 3, 1))
	r11 <- eiginv(qr.R(qrd)[1:p, 1:p])
	pe1 <- 0
	for (j in 1:p)
	{
		a[, , j] <- crossprod(r11, a[, , j]) %*% r11 * sigma
		pe1 <- pe1 + 2 * sum(a[, , j] ^ 2) + sum(diag(as.matrix(a[, , j]))) ^ 2
	}
	
	in1 <- 0
	for (j in (p + 1):rnk)
	{
		a[, , j] <- crossprod(r11, a[, , j]) %*% r11 * sigma
		in1 <- in1 + 2 * sum(a[, , j] ^ 2) + sum(diag(as.matrix(a[, , j]))) ^ 2
	}
	
	pe1 <- sqrt(pe1 / (p * (p + 2)))
	in1 <- sqrt(in1 / (p * (p + 2)))
	pe <- pe1 * sqrt(qf(19 / 20, p, n - p))
	int <- in1 * sqrt(qf(19 / 20, p, n - p))
	curvatures <- list(pe=pe, int=int, a=a,cutfp = 1/sqrt(qf(.95,p,n-p)))
	return(curvatures)
}
