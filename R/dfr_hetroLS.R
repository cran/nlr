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
dfr.hetroLS <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=1e-4, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),varmodel,tau=getInitial(varmodel,vdata),... ){

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	stage1 <- nlsnm(formula=formula,data=data,start=start,control=control)
	if(is.Fault(stage1)) return(stage1)
	t <- predict(stage1,newdata=data)
	ri <- residuals(stage1)
	n <- length(ri)
	nrp <- nonrepl(list(x=data$xr,y=data$yr))
	fdata <- NULL
	fdata[[formula$independent]] <- nrp$x[nrp$xm]
	fdata[[formula$dependent]] <- nrp$y[nrp$xm]
	z <- zvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]        # variance z=zi , si^2                #
	vdata<-as.list(NULL)
	vdata[[varmodel$dependent]] = z[nrp$xm]
	if(is.null(data[[varmodel$independent]])) 
		vdata[[varmodel$independent]] <- predict(stage1,newdata=fdata)                 # non replicated    #
	else
		vdata[[varmodel$independent]] <- data[[varmodel$independent]][nrp$xm]  # non replicated    #
	wi <- pmax(1,nrp$ni-1)
	start.tau <- tau
	data2 <-c(vdata,tau)
	varcomp <- eval(varmodel,data2)
	g <- ( 2.0 * as.numeric(varcomp$predictor)^2 / wi ) #/ stage1$parameters$sigma
	vm <- diag(g)
	rm <- diag(1.0/sqrt(g))

	#---------------------------------- step 2 ------------------
	stage2 <- nlsnm(formula=varmodel,data=vdata,start=tau,control=control,vm=vm,rm=rm)
	if(is.Fault(stage2)) return(stage2)

	ps <- stage1$parameters[names(formula$par)]
	vcn <- predict(stage2,newdata=vdata)
	vc <- rep(0,n)
	for (i in nrp$xo){
		vc[nrp$xo==i]=vcn[i]
	}
	g <- vc#/stage2$parameters$sigma^2
	vmat <- diag(g)
	sqg <- diag(sqrt(g))
	umat <- diag(sqg)
	tumat <- t(umat)
	rmat <- diag(1.0 / sqrt(g)) #eiginv(tumat,stp=F)
	if(is.Fault(rmat)) return(rmat)
	# ------------------------------- step 3 ----------------------------
	stage3 <- nlsnm(formula=formula,data=data,start=start,control=nlr.control(tolerance=tolerance*10, maxiter=maxiter,trace=trace,minlanda=minlanda),vm=vmat,rm=rmat)

	if(is.Fault(stage3)) result <- nl.fitt.gn(
						parameters =  stage1$parameters,
						correlation =  stage1$correlation ,
						form =         stage1$form ,
						response =     stage1$response,
						predictor =    stage1$predictor, 
						curvature =    stage1$curvature ,
						history =      stage1$history ,
						method = 		 stage1$method,
						data =         stage1$data ,
						sourcefnc =    match.call(),
						Fault =        stage1$Fault,
						vm=            vmat,
						rm=            rmat)

	else{
		result <- stage3
		result@method =	fittmethod(methodID=			1,
										methodBR=			14,
										detailBR=			"Modified newton",
										subroutine=		"nl.hetroLS")
	} 
	nvdata<-as.list(NULL)
	nvdata[[varmodel$independent]]<-vdata[[varmodel$independent]][nrp$xo]
	nvdata[[varmodel$dependent]]<-vdata[[varmodel$dependent]][nrp$xo]
	result@hetro<-nl.fitt(
								parameters=    stage2$parameters,
								form=          varmodel,
								predictor =    stage2$predictor, 
								response =     stage2$response, 
								history =      stage2$history, 
								method =       stage2$method,
								data =         nvdata,
								sourcefnc =    stage2$sourcefnc,
								Fault =        Fault(),
								#others =       list(refvar=stage2$objfnc$refvar)
								)

	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'dfr.hetro'                                 |
#|                                                                                 |
#|                           15   Sep 2013                                         |
#|                                                                                 |
#|                    Hossein Riazoshams, dep stat. stockholm U                    |
#|                                                                                 |
#+#################################################################################+
