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
nl.robhetro <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.000010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),robfunc,varmodel,tau=NULL,...){

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace
	FT <- NULL

	stage1 <- nlmest.NLM(formula, data=data, start=start,control=control,robfunc=robfunc,robscale=control$robscale,...)
	if(is.Fault(stage1)){
		print("nl.robhetro stage 1 first try failed, NLM is gone be evaluated.")
		stage1 <- nlmest.WF(formula, data=data, start=start,control=control,robfunc=robfunc,robscale=control$robscale,...)
		if(is.Fault(stage1)){
		  stage1 <- nlmest.NLM.sCase2(formula, data=data, start=start,control=control,robfunc=robfunc,robscale=control$robscale,...)    
		  if(is.Fault(stage1)){
        print("nl.robhetro stoped at stage 1 with error.")
			  return(stage1)
		  }
		}
	}
	if(trace) plot(stage1)
  
  ri <- residuals(stage1)
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	z <- rzvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]              ## variance z=zi , si^2          ##
	
	data2 <- NULL
	data2[[varmodel$dependent]] <- z[nrp$xm]              ##  vr, (nonreplicate)           ##
	if(is.null(data[[varmodel$independent]])) {
		t <- predict(stage1,newdata=stage1$data)
		data2[[varmodel$independent]] <- t[nrp$xm]        ##  t=mu     (nonreplicate)      ##
	}
	else
		data2[[varmodel$independent]] <- data[[varmodel$independent]]
	if(any(data2[[2]]<0)){
		return(Fault(FN=20))
	}


	###### stage 2 iterate  1    ########################################################################

	if(is.null(tau)) 
		if(is.null(varmodel$selfStart)) start.tau <- varmodel$par   # in R.. [[2]] is diferent from splus
		else{
			data2 <- NULL
			data2[[varmodel$dependent]] <- z[nrp$xm]
			if(is.null(data[[varmodel$independent]])) 
				data2[[varmodel$independent]] <- predict(stage1,newdata=stage1$data)[nrp$xm]
			start.tau <- getInitial(varmodel,data2)
		}
	else start.tau <- tau

#	print("start stage 2222--------------------------------------------------")

	stage2<- optim.NLM(objfnc=loss.robchis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,varmodel=varmodel,
		robfunc=robfunc,control=nlr.control(tolerance=tolerance,maxiter=maxiter*2,minlanda=minlanda),...)
	if(is.Fault(stage2)){#  || stage2$parameters$sg < 0){
			print("error in first try stage 2 occured, will try another starting in nl.robhetro")
			stage2<- optim.WF(objfnc=loss.robchis,data=data,start=start.tau,formula=formula,theta=stage1$parameters,varmodel=varmodel,
				robfunc=robfunc,control=nlr.control(tolerance=tolerance,maxiter=maxiter*2,minlanda=minlanda),...)
			if(is.Fault(stage2)){
				print("error at stage 2 robhetro")
				 return(stage2)
			}
	}
#	if(stage2$parameters$sg < 0){
#		print(stage2)
#		result <- stage2
#		result$Fault <- Fault(FN=19,FF="nl.robhetro")
#		return(result)
#	}

	###### stage 3    ########################################################################	
	ps <- stage1$parameters[names(formula$par)]
	vc <- as.numeric(stage2$objfnc$varcomp$predictor)
	gvar <- vc#/stage2$parameters$sg^2	
	vmat <- diag(gvar)
	rmat <- diag(1.0/sqrt(gvar))
#	ps["sigma"] <- stage2$parameters$sg
# print("stage 3-------------------=---------------------------------------------------------------------------------------")
	stage3 <- nlmest.NLM(formula, data=data, start=ps,
		control=nlr.control(tolerance=tolerance,maxiter=maxiter*10,minlanda=minlanda,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
	result <- stage3
	if(is.Faultwarn(stage3)){
		FT <- " error in stage 3 of nl.robhetro, trying another initial value"
		print("error in first try stage 3 occured, will try another starting in nl.robhetro")
		ps<-start
		ps["sigma"] <- stage2$parameters$sg
    stage3 <- nlmest.WF(formula, data=data, start=ps,
			control=nlr.control(tolerance=tolerance,maxiter=maxiter*3,minlanda=minlanda,robscale=control$robscale),vm=vmat,rm=rmat,robfunc=robfunc,...)
		if(is.Fault(stage3)) return(stage3)
		result <- stage3
	}
	htheta <- stage2$objfnc$value
	result@method<- fittmethod(methodID=5,subroutine ="nl.robhetro",lossfunction="robloss.gn")
	result@hetro<-nl.fitt.rob(
								parameters=    stage2$parameters,
								form=          varmodel,
								predictor =    stage2$objfnc$varcomp$predictor, 
								response =     stage2$objfnc$zi,
								history =      stage2$history, 
								method =       stage2$method,
								data =         stage2$objfnc$vcmdata,
								sourcefnc =    stage2$objfnc$sourcefnc,
								Fault =        Fault(FT=FT),
								htheta =       stage2$objfnc$value,
								rho =          stage2$objfnc$rho)
	result@others=list(refvar=stage2$objfnc$refvar)
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
