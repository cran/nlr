
#**********************************
#**  transform rmat * value
#**********************************
transform <- function(value,rm){
	g <- attr(value, "gradient")
	h <- attr(value, "hessian")
  value <- rm%*%value
  if(! is.null(g)) attr(value, "gradient") <- rm %*% g
	if(! is.null(h)) attr(value, "hessian") <- h %3d*m%rm
	return(value)
}
#**********************************
#**  transform inverse rmat(-1) * value
#**********************************
transforminv <- function(value,rm){
	rmatin <- eiginv(rm)

	return(transform(value,rm ))
}

#****************************************************************
#**  transquad calculate standard deviation and all derivatives
#** This is transforming from variance to standard deviation  
#****************************************************************
transfquad <- function(varcomp){
	vc <- as.numeric(varcomp)
	sd <- sqrt(vc)                                #  sigma i (standard deviation)      #
	if(! is.null(attr(varcomp,"gradient"))){
#				attr(sd,"gradient") <- (0.5*sd) * attr(varcomp,"gradient")         #  grad (sigmai)   gradient stdev    #
				attr(sd,"gradient") <- -attr(varcomp,"gradient")/(2*sd) 
	}
	if(! is.null(attr(varcomp,"hessian"))){
		vcgt <- t(attr(varcomp,"gradient"))
		gtg <- vcgt %m3d% attr(varcomp,"gradient")      # vcgT * vc (n*p*p)                  #
		.temp2 <- 4.0 * vc * sd                         # 4* sigma ^3                        #
		.expr2 <- (1.0 / .temp2) * gtg                  # 1/4* sigma^3 * gr T * grad         #
		.expr1 <- attr(varcomp,"hessian") / (2 * sd)    # 1/2*sigma  * hess                  #
		sdh <- .expr2 - .expr1
		attr(sd,"hessian") <- sdh
	}
	return(sd)
}

#****************************************************************
#**  transquadvec calculate standard deviation and all derivatives
#** This is transforming from variance to standard deviation
#** Inportant: diferent with transquad is: vector only not matrix.  
#****************************************************************
transfquadvec <- function(varcomp){
	vc <- as.numeric(varcomp)
	sd <- sqrt(vc)                                #  sigma i (standard deviation)      #
	if(! is.null(attr(varcomp,"gradient"))){
#				attr(sd,"gradient") <- (0.5*sd) * attr(varcomp,"gradient")         #  grad (sigmai)   gradient stdev    #
				attr(sd,"gradient") <- -attr(varcomp,"gradient")/(2*sd) 
	}
	if(! is.null(attr(varcomp,"hessian"))){
		vcgt <- attr(varcomp,"gradient") ^ 2
		.temp2 <- 4.0 * vc * sd                         # 4* sigma ^3                        #
		.expr2 <- (1.0 / .temp2) * vcgt                 # 1/4* sigma^3 * gr T * grad         #
		.expr1 <- attr(varcomp,"hessian") / (2 * sd)    # 1/2*sigma  * hess                  #
		sdh <- .expr2 - .expr1
		attr(sd,"hessian") <- sdh
	}
	return(sd)
}
