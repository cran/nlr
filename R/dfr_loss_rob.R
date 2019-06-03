#************************************************************************
#**   robust loss unction                                              **
#**     start: argument must be a list containing "sigma".             **
#**     argument ... is the arguments for 'robfunc' extra argument     **
#**     for 'robfunc' the robust rho function.                         **
#**     Fault: = 9,  sigma is not a part of start list.                **
#************************************************************************
dfr.robloss <- function(formula,data,start,robfunc,control=nlr.control(),rmat=NULL,...){
	Fault2 <- Fault()
	data <- as.list(data)
	start <- as.list(start)
	.tmp <-  (is.null(start$sigma))
	if(.tmp) return(list(Fault=Fault(FL=T,FN=8,FF="robloss")))
	datalist <- c(start,data)
	fmod <- eval(formula,datalist)						#  f(t) n*1									#
	if(is.Fault(fmod)) return(fmod)
	####################### GN
	if(! is.null(rmat)) {
		fmod$predictor <- transformNR(fmod$predictor,rmat)
		fmod$response <- transformNR(fmod$response,rmat)
	}
	###########################	
	prd <- as.numeric(fmod$predictor)					#  f(t) values								#
	n <- length(prd)
	p <- formula$p
	gprd <- attr(fmod$predictor,"gradient")			#  f'(t)	Gradient	n*p						#
	hprd <- attr(fmod$predictor,"hessian")			#  f"(t)	Hessian	n*p*p					#
	rsd <- as.numeric(fmod$response) - prd			#  ri		Residals	n*1						#
	gres <- -1 * gprd                           #  J=r'i=-grad(f)		n*p						#
	hres <- -hprd                               #  r"i					n*p*p					#
	attr(rsd,"gradient") <- gres
	attr(rsd,"hessian") <- hres	
  ris <- rsd / start$sigma 
	if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
	else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
	rho <- as.numeric(robvalue)
	htheta <- sum(rho)

	## ----------------------------------------- ##
	##     numerical derivative                  ##
	## ----------------------------------------- ##

	xnames <- names(start)
	x <- unlist(start[names(start) != "sigma"])
	x0 <- x
   fvec <- htheta                         # f(x), if is NULL compute here.
	eps <- sqrt(.Machine$double.eps)
	fgrad <- matrix(rep(0,p),nrow=p)                                ## fgrad is gradient (jacobian) #
	fhessian <- matrix(rep(0,p*p),nrow=p)
	for(j in 1:p){
		temp <- x[j]
		h <- eps*abs(temp)
		if(h == 0) h <- eps
		x[j] <- temp + h
		t <- as.list(x)                                              ## calculate function at xj+h
		t$sigma <- start$sigma
		datalistx <- c(data,t)
		objw <- eval(formula,datalistx)
		rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
		if(! is.null(rmat)) rix <- rmat %*% rix
		ris <- rix / start$sigma 
		if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
		else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
		rho <- as.numeric(robvalue)
		wa <- sum(rho)

		x[j] <- temp - h                                             ## centered operator
		t <- as.list(x)                                              ## calculate function at xj+h
		t$sigma <- start$sigma
		datalistx <- c(data,t)
		objw <- eval(formula,datalistx)
		rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
		if(! is.null(rmat)) rix <- rmat %*% rix
		ris <- rix / start$sigma 
		if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
		else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
		rho <- as.numeric(robvalue)
		wa2 <- sum(rho)
		fgrad[j] = (wa - wa2)/ (2*h)                                 ## centered          #
		x[j] <- temp - 2*h                                           ## x-2h
		t <- as.list(x)                                              ## calculate function at xj+h
		t$sigma <- start$sigma
		datalistx <- c(data,t)
		objw <- eval(formula,datalistx)
		rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
		if(! is.null(rmat)) rix <- rmat %*% rix
		ris <- rix / start$sigma 
		if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
		else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
		rho <- as.numeric(robvalue)
		wa3 <- sum(rho)

		x[j] <- temp + 2*h                     ### x+2h
		t <- as.list(x)                                              ## calculate function at xj+h
		t$sigma <- start$sigma
		datalistx <- c(data,t)
		objw <- eval(formula,datalistx)
		rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
		if(! is.null(rmat)) rix <- rmat %*% rix
		ris <- rix / start$sigma 
		if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
		else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
		rho <- as.numeric(robvalue)
		wa4 <- sum(rho)

		fhessian[j,j] <- (-wa4 + 16*wa-30*fvec+16*wa2-wa3) / (12*h*h)
		i <- j+1
		while(i <= p){
			temp2 <- x[i]
			h2 <- eps*abs(temp2)
			if(h2 == 0) h2 <- eps

			x[i] <- temp2 + h2
			x[j] <- temp + h                     # x+hiei+hjej   in SAS notation not our here
			                                     # important note i,j switch here
			t <- as.list(x)                                              ## calculate function at xj+h
			t$sigma <- start$sigma
			datalistx <- c(data,t)
			objw <- eval(formula,datalistx)
			rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
			if(! is.null(rmat)) rix <- rmat %*% rix
			ris <- rix / start$sigma 
			if(control$robscale) robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
			else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
			rho <- as.numeric(robvalue)
			wb1 <- sum(rho)

			x[i] <- temp2 - h2
			x[j] <- temp + h                     # x+hiei-hjej
			t <- as.list(x)                                              ## calculate function at xj+h
			t$sigma <- start$sigma
			datalistx <- c(data,t)
			objw <- eval(formula,datalistx)
			rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
			if(! is.null(rmat)) rix <- rmat %*% rix
			ris <- rix / start$sigma 
			if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
			else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
			rho <- as.numeric(robvalue)
			wb2 <- sum(rho)

			x[i] <- temp2 + h2
			x[j] <- temp - h                     # x-hiei+hjej
			t <- as.list(x)                                              ## calculate function at xj+h
			t$sigma <- start$sigma
			datalistx <- c(data,t)
			objw <- eval(formula,datalistx)
			rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
			if(! is.null(rmat)) rix <- rmat %*% rix
			ris <- rix / start$sigma 
			if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
			else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
			rho <- as.numeric(robvalue)
			wb3 <- sum(rho)

			x[i] <- temp2 - h2
			x[j] <- temp - h                     # x-hiei-hjej
			t <- as.list(x)                                              ## calculate function at xj+h
			t$sigma <- start$sigma
			datalistx <- c(data,t)
			objw <- eval(formula,datalistx)
			rix <- as.numeric(objw$response) - as.numeric(objw$predictor)
			if(! is.null(rmat)) rix <- rmat %*% rix
			ris <- rix / start$sigma 
			if(control$robscale)  robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)
			else robvalue <- robfunc$fnc(ris, ...)			#  rho(ri/sigma)	n*1							#
			rho <- as.numeric(robvalue)
			wb4 <- sum(rho)
			numerator <- (wb1-wb2-wb3+wb4)
			denomerator <- (4*h*h2)

			fhessian[i,j] <- fhessian[j,i] <- numerator /denomerator     ### centered
			i <- i+1
		}
		x <- x0
	}
	attr(htheta,"gradient") <- fgrad
	attr(htheta,"hessian") <- fhessian
	## ----------------------------------------- ##
	##   End of  numerical derivative            ##
	## ----------------------------------------- ##

	result <- list(htheta=htheta,rho=robvalue,ri=rsd,fmod=fmod,Fault=Fault2)
	#****************************************************
	#*****   output:                                ***** 
	#*****  htheta=sum(rho)  fixed number           *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  rho(ri/sigma)       n*1 vector          *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  ri               residuals              *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  h(theta)=he	ss.... p1+p2 hessian of (h)*****

	#*****  fmod:  computed function contains       *****
	#*****  response and or its gradient and hessian*****
	#*****  predictor and or its gradient & hessian *****
	#****************************************************

	return(result)
}
