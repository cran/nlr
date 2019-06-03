
#+------------------------------------------------------------------------------+
#|    nloutJL: find the measurments for outliers, using Jacobian Leverage.      |
#|       nlfited: nl.fitt or nl.fitt.rob or nl.fitt.mm or any inherited object. |
#|    detail. the residuals and x can be robust but it make confusion.          |
#|    this function is only asmall modification of nlout.                       |
#+------------------------------------------------------------------------------+
nlout.JL<-function(nlfited){
	if(class(nlfited)=="Fault") return(nlfited)
	if(nlfited$Fault$FN==T) return(nlfited$Fault)
	v <- attr(nlfited@predictor,"gradient")
	hs <- attr(nlfited@predictor,"hessian")
	res <- residuals(nlfited)
	p <- ncol(v)
	n <- length(res)
	dataset <- nlfited$data[c(nlfited$form$independent,nlfited$form$dependent)]
	dataset <- data.frame(dataset)
	JL <- jaclev(v,hs,res)
	if(class(JL)=="Fault") return()
	jvmat <- diag(JL)
  sigma <- nlfited$scale
	jl.studres <- res / sigma / sqrt(1-jvmat)				##  studentized residual, ri / (sigm sqrt(1-wii))
	jl.elliptnorm <- jl.studres^2 * jvmat / (1-jvmat) / p	##  Influence curve, Elliptic norm (Cook d)
	jl.hadi <- jvmat / (1-jvmat)								##  Hadi to assess high leverage points wii/(1-wii)
	mdhad <- median(jl.hadi)
	jl.ctfhadi1 <- mdhad + 2 * mad(jl.hadi)
	jl.ctfhadi2 <- mdhad + 3 * mad(jl.hadi)

	###############################################
   ##      Single deletion                      ##
	###############################################
	d.stud <- d.sigma <- d.yhat <- d.ffits <- rep(0,n)    
	yi <- as.numeric(nlfited$predictor)
	yihat <- as.numeric(nlfited$response)
			############################################### new idea for deleted object   #############
	for(i in 1:n){
		dst2 <- dataset[-i,]									##   xy(-i)
		newfiti <- recalc(nlfited,dst2)					##   fit all (-i)
		newfiti@Fault@FF <- "nlout"
		if(newfiti$Fault$FL) return(newfiti$Fault)
		d.sigma[i]  <- newfiti$parameters[["sigma"]] 	##  sigma(-i)
		yhati <- predict(newfiti,newdata=dataset[i,])				##  yhat(-i) grd&hess
		d.yhat[i] <- as.numeric(yhati)						##  yhat(-i)
	}
	jl.delstud <- (yi - d.yhat) / d.sigma / sqrt(1-jvmat)	##  deleted resid
	yhdiff <- yihat - d.yhat								##  yi^ - y(-i)^
	jl.d.ffits <- yhdiff / d.sigma / sqrt(jvmat) 				##  dffits
	ctfdffits <- 2 * sqrt(p/n)								##  cut of poin dffits
	jl.atkinson <- sqrt( (n-p) * jvmat / p / (1-jvmat) ) * abs(jl.d.ffits)   ##  Atkinsons distance
	################################################################### output   #########################
	result1 <- structure(.Data=list(
								jvmat,
								nl.robmeas(measure=jl.studres,cutofpoint=c(2.5,3,-2.5,-3),name="Jacobian Leverage Studentised Residuals"),
								nl.robmeas(measure=sqrt(jl.elliptnorm),cutofpoint=1,name="Jacobian Leverage Elliptic Norm (Cook Dist)"),
								nl.robmeas(measure=jl.hadi,cutofpoint=c(jl.ctfhadi1,jl.ctfhadi2),name="Jacobian Leverage Hadi potential"),
								nl.robmeas(measure=jl.delstud ,cutofpoint=c(2.5,3),name="Jacobian Leverage Deletion Studentized"),
								nl.robmeas(measure=jl.d.ffits ,cutofpoint=ctfdffits ,name="Jacobian Leverage DFFITS"),
								nl.robmeas(measure=jl.atkinson ,cutofpoint=2,name="Jacobian Leverage Atkinson Distance")
							  ),
							.Names=c(
										"jl.vmat", 
										"jl.studres",
										"jl.cook",
										"jl.hadi",
										"jl.delstud",
										"jl.dffits",
										"jl.atk")
							)
	return(result1)
}

#--------------------------------------------------------------------------------
#|                                  End of nlout. 
#|                                  31/12/2008
#--------------------------------------------------------------------------------

