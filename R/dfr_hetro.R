
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'dfr.hetro', cllasic estimate of a nonlinear function.      | *
#* |       with hetro variance model function. derivative free.             | *
#* |   generalized                                                          | *
#* |  Note: becarefull to using this function when there is not outlier, it | *
#* |    may not work witout outlier, in this case better to use nlmest      | *
#* |  the problem is in part of two (p2) in hessian its big here.           | *
#* |    argumnts:                                                           | *
#* |      formula: 'nl.form' object, the function mode.                     | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      start:   starting values, it must contains 'sigma', selstart      | *
#* |               for nl.form object is not created yet, take cre of it.   | *
#* |      varfnc: var function, it must be nl.form of variance models       | *
#* |      tau:     starting value of tau. if is null the stored value in    | *
#* |               vardnc object of nl.form will be stored.                 | *
#* |      ...:     can be entries for robust loss function parameters.      | *
#* |                                                                        | *
#* |   Important Note: variance must be a product function in sigma, i.e.   | *
#* |       varfunc = sigma^2 * h(f,tau)                                     | *
#* |       in feature the general form will be added.                       | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
dfr.hetro <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.000010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),varmodel,tau=NULL,... ){

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	#on.exit(print(stage2$parameters),frame=1)
	FT <- NULL

	stage1 <- nlsnm(formula=formula,data=data,start=start,control=nlr.control(tolerance=tolerance,maxiter=10*maxiter),...)
	if(is.Fault(stage1)){
		print ("stage 1 hetro")
		return(stage1)
	}
	ri <- residuals(stage1)
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	z <- zvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]              # variance z=zi , si^2          #
	data2 <- NULL

	#*********************** stage 2.....................
	
	if(is.null(tau)) 
		if(is.null(varmodel$selfStart)) start.tau <- varmodel$par  # different splus
		else{
			data2 <- NULL
			data2[[varmodel$dependent]] <- z[nrp$xm]
			if(is.null(data[[varmodel$independent]])) 
				data2[[varmodel$independent]] <- as.numeric(predict(stage1,newdata=stage1$data))[nrp$xm]
			start.tau <- getInitial(varmodel,data2)
		}
	else start.tau <- tau
	stage2 <- optim.NM(objfnc=loss.chis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,
		varmodel=varmodel,tolerance=tolerance,...)

	if(is.Fault(stage2)){ # || stage2$parameters$sg < 0){
		 return(stage2)
	}

	# ***************** stage 3 .....................

	ps <- stage1$parameters
	vc <- as.numeric(stage2$objfnc$varcomp$predictor)
	gvar <- vc#/stage2$parameters$sg^2
	vmat <- diag(gvar)
	rmat <- diag(1.0/sqrt(gvar))

	stage3 <- nlsnm(formula, data=data, start=ps,tolerance=tolerance, 
		maxiter=maxiter,vm=vmat,rm=rmat)

	if(is.Faultwarn(stage3)){
		ps <- start
		return(stage3)
	}
	result <- stage3
	result@method = 		 fittmethod(methodID=			4,
											methodBR=			13,
											subroutine=		"dfr.hetro")
	result@hetro<-nl.fitt(
								parameters=    stage2$parameters,
								form=          varmodel,
								predictor =    stage2$objfnc$varcomp$predictor, 
								response =     stage2$objfnc$zi, 
								history =      stage2$history, 
								method = 		 fittmethod(methodID=			4,
																methodBR=			13,
																subroutine=		"optim.NM"),
								data =         stage2$objfnc$vcmdata,
								sourcefnc =    stage2$objfnc$sourcefnc,
								Fault =        Fault(FT=FT),
								others =       list(refvar=stage2$objfnc$refvar))
	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'dfr.hetro'                                 |
#|                                                                                 |
#|                              Aug 2013                                           |
#|                                                                                 |
#|                    Hossein Riazoshams, Stockholm U                              |
#|                                                                                 |
#+#################################################################################+
