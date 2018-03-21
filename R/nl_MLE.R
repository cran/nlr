#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nl.MLE', MLE estimate of a                                 | *
#* |       nonlinear function. with hetro variance model function,          | *
#* |       and weights.                                                     | *
#* |                                                                        | *
#* |    Exactly same as WM-est, but for fast computing writen again.        | *
#* |    argumnts:                                                           | *
#* |      formula: 'nl.form' object, the function mode.                     | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      start:   starting values, it must contains theta.                 | *
#* |               for nl.form object is not created yet, take cre of it.   | *
#* |      vm,rm,   weights, this have priority, if provided varmodel will   | *
#* |               be ignored.                                              | *
#* |      varmodel: var function, it must be nl.form of variance models     | *
#* |      tau:     starting value of tau. if is null the stored value in    | *
#* |               vardnc object of nl.form will be stored.                 | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
nl.MLE <-  function(formula, data, start=getInitial(formula,data),vm=NULL,rm=solve(t(chol(vm))),control=nlr.control(derivfree=T)
		,varmodel=NULL,tau=varmodel$par,...){

	if(! is.null(tau)) startlist <- c(start,tau)
	derivfree <- control$derivfree
  startlist <- as.list(startlist)

	if (is.null(startlist$sigma)){
		dt1 <- c(data,startlist)
		ht <- eval(formula,dt1)
		dv <- as.numeric(ht$predictor) - as.numeric(ht$response)
		if(! is.null(vm)) dv <- rm %*% dv		
		startlist[["sigma"]] <- sd(dv)   # mscale(dv)
	}
	# +----------------------------------------------------------------+ *
	# |     derivative free                                            | *
	# +----------------------------------------------------------------+ *
	if(derivfree){
		if(is.null(vm))
						if(is.null(varmodel)) result <- nlsnm(formula=formula,start=start,control=control,...)
				else{                             #********** mle der free, Nelder Mead *****************
            control$derivfree=T           #........... only now later derivfree must be removed
						result <- nl.robhetroWM(formula=formula,data=data,start=startlist,varmodel=varmodel,
							robfunc=nlr::nl.robfuncs[["least square"]],control=control,tau=tau,...)
				}
		else
				result <- nlsnm(formula=formula, data=data,start=start,control=control,vm=vm,rm=rm,...)
	}
	# +----------------------------------------------------------------+ *
	# |     derivative exist                                           | *
	# +----------------------------------------------------------------+ *
	else{		
		if(is.null(vm))
					if(is.null(varmodel)) 
						result <- nlsqr(formula=formula,data=data,start = start,control=control,...)
				else
						result <- nl.robhetroWM(formula=formula,data=data,start=startlist,varmodel=varmodel,
							robfunc=nlr::nl.robfuncs[["least square"]],control=control,tau=tau,...)
		else
				result <- nlsqr.gn(formula=formula,data=data,start = start,control=control,vm=vm,rm=rm,...)	
	}
	result@method@methodID <- 12
	result@method@method <- "MLE"
	result@method@detail <- "Maximum Likelihood Estimate"
	
	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.MLE                                     |
#|                                                                                 |
#|                             22 Aug 2013                                         |
#|                                                                                 |
#|                    Hossein Riazoshams, Stat Dep SU                              |
#|                                                                                 |
#+#################################################################################+



