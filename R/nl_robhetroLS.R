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
#* |      Method:  is to control when error happens and need control        | *
#* |               manually, program try to control errors but other        | *
#* |               un predicted error like log(0) may happens.              | *
#* |                                                                        | *
#* |   Important Note: variance must be a product function in sigma, i.e.   | *
#* |       varfunc = sigma^2 * h(f,tau)                                     | *
#* |       in feature the general form will be added.                       | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
nl.robhetroLS <-  function(formula, data, start=getInitial(formula,data),
	control=nlr.control(tolerance=0.00001, minlanda=1 / 2 ^ 10, maxiter=30 * length(start),robscale=T),robfunc,varmodel,
	tau=varmodel$par,...){
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	stage1 <- nlmest.NLM(formula, data=data, start=start,robfunc=robfunc,control=control,...)
	if(is.Fault(stage1)) return(stage1)
	if(control$trace) print(stage1$parameters)
	t <- predict(stage1,newdata=data)
	ri <- residuals(stage1)
	n <- length(ri)
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	fdata <- NULL
	fdata[[formula$independent]] <- nrp$x[nrp$xm]
	fdata[[formula$dependent]] <- nrp$y[nrp$xm]
	z <- rzvalues(ri,nrp$ni,nrp$xo)#[nrp$xm]        # variance z=zi , si^2                #
	vdata<-as.list(NULL)
	vdata[[varmodel$dependent]] = z[nrp$xm]
	if(is.null(data[[varmodel$independent]])) 
		vdata[[varmodel$independent]] <- predict(stage1,newdata=fdata)                      # non replicated      #
	else
		vdata[[varmodel$independent]] <- data[[varmodel$independent]][nrp$xm]       # non replicated      #
	wi <- pmax(1,nrp$ni-1)
	vm <- diag(wi)
	rm <- diag(1.0/sqrt(wi))
	###################################### step 1	iteration 1  ############
	startv <- tau
	data2 <-c(vdata,tau)
	varcomp <- eval(varmodel,data2)
	g <- 2.0 * as.numeric(varcomp$predictor)^2 / wi / tau$sg^2   ### then v(Rzi)=sg^2
	vm <- diag(g)
	rm <- diag(1.0/sqrt(g))
#	print("start stage 2---1--===============================================================")
	stage2 <- nlmest.WF(varmodel, data=vdata, start=tau,control=nlr.control(tolerance=tolerance, minlanda=minlanda, maxiter=5*maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
	if(is.Faultwarn(stage2) ||  stage2$parameters$sg<0){
#		print("error in first try stage2..1... will try new initial value nl.robhetroLS")
		start.tau <- getInitial(varmodel,vdata)
		stage2 <- nlmest.NLM(varmodel, data=vdata, start=tau,control=nlr.control(tolerance=tolerance, minlanda=minlanda, maxiter=5*maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
		if(is.Fault(stage2)) return(stage2)
	}
	if(control$trace) print(stage2$parameters)
	###################################### step 2    iteration 2#############
#	print("second iteration step 2-2...............................................")
	startv <- stage2$parameters[names(varmodel$par)]
	data2 <-c(vdata,startv)
	varcomp <- eval(varmodel,data2)
	g <- 2.0 * as.numeric(varcomp$predictor)^2 / wi / startv$sg^2
	vm <- diag(g)
	rm <- diag(1.0/sqrt(g))
	stage2 <- nlmest.WF(varmodel, data=vdata, start=startv,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=5*maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
	if(is.Fault(stage2) ||  stage2$parameters$sg<0){
		print("error in first try stage2..2... will try new initial value nl.robhetroLS")
		stage2 <- nlmest.NLM(varmodel, data=vdata, start=startv,control=nlr.control(tolerance=tolerance*5, minlanda=minlanda, maxiter=5*maxiter,trace=trace,robscale=control$robscale),robfunc=robfunc,vm=vm,rm=rm,...)
		if(is.Fault(stage2)) return(stage2)
	}
	if(control$trace) print(stage2$parameters)
	############################################ stage 3
	ps <- stage1$parameters[names(formula$par)]
	vcn <- predict(stage2,newdata=vdata)
	vc <- rep(0,n)
	for (i in nrp$xo){
		vc[nrp$xo==i]=vcn[i]
	}
	g <- (vc/stage2$parameters$sg^2)
	vmat <- diag(g)
	umat <- diag(sqrt(g))
	tumat <- t(umat)
	rmat <- diag(1.0 / diag(tumat)) #eiginv(tumat,stp=F)
	if(is.Fault(rmat)) return(rmat)
	stage3 <- nlmest.NLM(formula, data=data, start=ps,control=nlr.control(maxiter=5*maxiter,tolerance=tolerance*10,minlanda=minlanda/10,trace=trace,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
	if(is.Faultwarn(stage3)){
#		if(! is.Fault(stage3)) #ps<-stage3$parameters[names(formula$par)]
		print("Warning: stage 3 not converged, will try another method WF")
		stage32 <- nlmest.WF(formula, data=data, start=ps,control=nlr.control(tolerance=tolerance*10,minlanda=minlanda/10, maxiter=maxiter*5,trace=trace,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
			if(is.Faultwarn(stage32)){
				if(! is.Fault(stage32)) ps<-stage32$parameters[names(formula$par)]		
				else ps <- start
				stage3 <- nlmest.NLM(formula, data=data, start=start,
					control=nlr.control(tolerance=tolerance*20, minlanda=minlanda, maxiter=maxiter,trace=trace,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
				if(is.Fault(stage3) & (! is.Fault(stage32))) stage3<-stage32
			}
			else stage3 <- stage32
	}
	result <- stage3
	if(is.Fault(stage3)) return(stage3)
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
