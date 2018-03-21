#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nl.rhetro', MM robust estimate of a nonlinear function.    | *
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
#* |      varmodel: var function, it must be nl.form of variance models     | *
#* |      tau:     starting value of tau. if is null the stored value in    | *
#* |               vardnc object of nl.form will be stored.                 | *
#* |      ...:     can be entries for robust loss function parameters.      | *
#* |                                                                        | *
#* |   Important Note: variance must be a product function in sigma, i.e.   | *
#* |       varfunc = sigma^2 * h(f,tau)                                     | *
#* |       in feature the general form will be added.                       | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
dfr.robhetroLS <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.001, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),
    robfunc,varmodel,tau=varmodel$par,method="NM",...){
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace
	switch(method,
		"NLM"={
			stage1 <- dfrmest.NLM(formula, data=data, start=start,control=control,robfunc=robfunc,...)
			if(is.Fault(stage1)){
				print("nl.robhetro stoped at stage 1 with error.")
				stage1 <- nlmest.NM(formula, data=data, start=start,control=control,robfunc=robfunc,...)
				if(is.Fault(stage1)) return(stage1)
			}
		},
		"NM"={
			stage1 <- nlmest.NM(formula, data=data, start=start,control=control,robfunc=robfunc,...)
			if(is.Fault(stage1)) return(stage1)
		}
	)
	t <- predict(stage1,newdata=data)
	ri <- residuals(stage1)
	n <- length(ri)
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	fdata <- NULL
	fdata[[formula$independent]] <- nrp$x[nrp$xm]
	fdata[[formula$dependent]] <- nrp$y[nrp$xm]
	z <- rzvalues(ri,nrp$ni,nrp$xo)                                                # variance z=zi , si^2                #
	vdata<-as.list(NULL)
	vdata[[varmodel$dependent]] = z[nrp$xm]
	if(is.null(data[[varmodel$independent]])) 
		vdata[[varmodel$independent]] <- as.numeric(predict(stage1,newdata=fdata))          # non replicated      #
	else
		vdata[[varmodel$independent]] <- data[[varmodel$independent]][nrp$xm]       # non replicated      #
	wi <- pmax(1,nrp$ni-1)
	vm <- diag(wi)
	rm <- diag(1.0/sqrt(wi))
	###################################### stage 2 iterate 1
	startv <- tau
	data2 <-c(vdata,tau)
	varcomp <- eval(varmodel,data2)
	g <- 2.0 * as.numeric(varcomp$predictor)^2 / wi / stage1$parameters$sigma    ### then v(Rzi)=sg^2
	vm <- diag(g)
	rm <- diag(1.0/sqrt(g))
	switch(method,
		"NLM"={
			stage2 <- dfrmest.NLM(varmodel, data=vdata, start=tau,tolerance=tolerance, minlanda=minlanda, 
				axiter=maxiter,robfunc=robfunc,vm=vm,rm=rm,control=control,...)
			if(is.Fault(stage2)){
				print("nl.robhetro stoped at stage 2 with error.")
				stage2 <- nlmest.NM(varmodel, data=vdata, start=tau,control=control,robfunc=robfunc,vm=vm,rm=rm,...)
				if(is.Fault(stage2)) return(stage2)
			}
		},
		"NM"={
			stage2 <- nlmest.NM(varmodel, data=vdata, start=tau,control=control,robfunc=robfunc,vm=vm,rm=rm,...)
			if(is.Fault(stage2)) return(stage2)
		}
	)

	###################################### stage 2  iterate 2
	startv <- stage2$parameters[names(varmodel$par)]

	data2 <-c(vdata,startv)
	varcomp <- eval(varmodel,data2)
	g <- 2.0 * as.numeric(varcomp$predictor)^2 / wi / stage1$parameters$sigma
	vm <- diag(g)
	rm <- diag(1.0/sqrt(g))
	switch(method,
		"NLM"={
			stage2 <- dfrmest.NLM(varmodel, data=vdata, start=startv,tolerance=tolerance, minlanda=minlanda, 
				axiter=maxiter,robfunc=robfunc,vm=vm,rm=rm,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=maxiter,robscale=control$robscale),...)
			if(is.Fault(stage2)){
				print("nl.robhetro stoped at stage 2 with error.")
				stage2 <- nlmest.NM(varmodel, data=vdata, start=startv,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
				if(is.Fault(stage2)) return(stage2)
			}
		},
		"NM"={
			stage2 <- nlmest.NM(varmodel, data=vdata, start=startv,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
			if(is.Fault(stage2)) return(stage2)
		}
	)
	############################################ stage 3
	ps <- stage1$parameters[names(formula$par)]
	vcn <- predict(stage2,newdata=vdata)
	vc <- rep(0,n)
	for (i in nrp$xo){
		vc[nrp$xo==i]=vcn[i]
	}
	g <- vc / stage1$parameters$sigma^2 #(vc/stage2$parameters$sg^2)
	vmat <- diag(g)
	umat <- diag(sqrt(g))
	tumat <- t(umat)
	rmat <- diag(1.0 / diag(tumat)) #eiginv(tumat,stp=F)

	switch(method,
		"NLM"={
			stage3 <- dfrmest.NLM(formula, data=data, start=ps,
				control=nlr.control(maxiter=3*maxiter,minlanda=minlanda, tolerance=tolerance,trace=trace,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
			if(is.Fault(stage3)){
				print("nl.robhetro stoped at stage 1 with error.")
				stage3 <- nlmest.NM(formula, data=data, start=ps,control=nlr.control(tolerance=tolerance, minlanda=minlanda, maxiter=3*maxiter,trace=trace),vm=vmat,rm=rmat,robfunc=robfunc,...)
				if(is.Fault(stage1)) return(stage3)
			}
		},
		"NM"={
			stage3 <- nlmest.NM(formula, data=data, start=ps,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=3*maxiter,trace=trace),vm=vmat,rm=rmat,robfunc=robfunc,...)
			if(is.Fault(stage3)) return(stage1)
		}
	)


	result <- stage3
	result@method@methodID <- 2
	result@method@subroutine <- "nl.robhetroLS"
	nvdata<-as.list(NULL)
	nvdata[[varmodel$independent]]<-vdata[[varmodel$independent]][nrp$xo]
	nvdata[[varmodel$dependent]]<-vdata[[varmodel$dependent]][nrp$xo]
	result@hetro<-nl.fitt.rob(
								parameters=    stage2$parameters,
								form=          varmodel,
								predictor =    stage2$predictor, 
								response =     stage2$response, 								
								history =      stage2$history, 
								method = 		 stage2$method,
								data =         nvdata,
								sourcefnc =    stage2$sourcefnc,
								Fault =        Fault(),
								others =       list(refvar=stage2$others$refvar),
								htheta =       stage2$htheta,
								rho =          stage2$rho)
	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.robhetro'                               |
#|                                                                                 |
#|                              Nov 2009                                           |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
