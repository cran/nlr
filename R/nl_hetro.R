
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nl.hetro', cllasic estimate of a nonlinear function.       | *
#* |       with hetro variance model function.                              | *
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
nl.hetro <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.000010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),varmodel,tau=NULL,... ){
	#on.exit(print(stage2$parameters),frame=1)
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	FT <- NULL
	stage1 <- nlsqr(formula, data=data, start=start, control=control)

  if(is.Fault(stage1)){
		print ("stage 1 hetro is fault")
		stage1 <- nlsnm(formula=formula,data=data,start=start,control=nlr.control(tolerance=tolerance,maxiter=10*maxiter,minlanda=minlanda),...)
		if(is.Fault(stage1)) return(stage1)
	}
	ri <- residuals(stage1)
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	z <- zvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]              # variance z=zi , si^2          #
	data2 <- NULL

	#*********************** stage 2.....................
	if(is.null(tau)) 
		if(is.null(varmodel$selfStart)) start.tau <- varmodel$par    # in R only [[2]] is different with splus
		else{
			data2 <- NULL
			data2[[varmodel$dependent]] <- z[nrp$xm]
			if(is.null(data[[varmodel$independent]])) 
				data2[[varmodel$independent]] <- predict(stage1,newdata=stage1$data)[nrp$xm]
      start.tau <- getInitial(varmodel,data2)
		}
	else start.tau <- tau

	stage2 <- optim.NLM(objfnc=loss.chis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,
				varmodel=varmodel,control=control,...)
	if(is.Fault(stage2)){ # || stage2$parameters$sg < 0){
		print("in stage 2 of nl_hetro first try failed with start.tau")
		stage2 <- optim.WF(objfnc=loss.chis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,
		varmodel=varmodel,control=control,...)
		if(is.Fault(stage2)){
			stage2 <- optim.NM(objfnc=loss.chis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,
				varmodel=varmodel,control=control,...)
			
			print("error stage 2 hetro")
			if(is.Fault(stage2)) return(stage2)
		}
	}
#	if(stage2$parameters$sg < 0) return(Fault(FN=19))
		
	# ***************** stage 3 .....................

	ps <- stage1$parameters
#	ps$sigma <- stage2$parameters$sg
	vc <- as.numeric(stage2$objfnc$varcomp$predictor)
	gvar <- vc#/stage2$parameters$sg^2
	vmat <- diag(gvar)
	rmat <- diag(1.0/sqrt(gvar))
  stage3 <- nlsqr.gn(formula, data=data, start=ps,control=nlr.control(minlanda=minlanda,maxiter=maxiter*3,tolerance=tolerance)
		,vm=vmat,rm=rmat)

	if(is.Faultwarn(stage3)){
		ps <- start
#		ps["sigma"] <- stage2$parameters$sg
		stage3 <- nlsqr.gn(formula, data=data, start=ps,control=nlr.control(minlanda=minlanda,maxiter=maxiter*20,tolerance=tolerance)
			,vm=vmat,rm=rmat)
		FT <- "in nl.hetro stage 3 first try failed will try with another imnitial value"
		if(is.Fault(stage3)){
			stage3 <- nlsnm(formula, data=data, start=ps,control=control,vm=vmat,rm=rmat)
		}
		if(is.Fault(stage3)) return(stage3)
	}
	result <- stage3
	result@method = 		 fittmethod(methodID=		4,
							methodBR=		14,
							detailBR=		"Modified Newton",
							subroutine=		"nl.hetro")
	result@hetro<-nl.fitt(
								parameters=    stage2$parameters,
								form=          varmodel,
								predictor =    stage2$objfnc$varcomp$predictor, 
								response =     stage2$objfnc$zi, 
								history =      stage2$history, 
								method = 		 fittmethod(methodID=		4,
													methodBR=		14,
													detailBR=		"Modified Newton",
													subroutine=		"optim.NLM"),

								data =         stage2$objfnc$vcmdata,
								sourcefnc =    stage2$objfnc$sourcefnc,
								Fault =        Fault(FT=FT),
								others =       list(refvar=stage2$objfnc$refvar))
	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.hetro'                                  |
#|                                                                                 |
#|                              Sep 2009                                           |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
