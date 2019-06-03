#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nl.robhetroWM', Weighted M robust estimate of a            | *
#* |       nonlinear function. with hetro variance model function.          | *
#* |       Proposed by Lim                                                  | *
#* |    argumnts:                                                           | *
#* |      formula: 'nl.form' object, the function mode.                     | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      start:   starting values, it must contains theta.                 | *
#* |               for nl.form object is not created yet, take cre of it.   | *
#* |      varmodel: var function, it must be nl.form of variance models     | *
#* |      tau:     starting value of tau. if is null the stored value in    | *
#* |               vardnc object of nl.form will be stored.                 | *
#* |      delta    optional entry can be used for segmenting initials       | *
#* |      ...:     can be entries for robust loss function parameters.      | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
nl.robhetroWM <-  function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.0001, minlanda=1 / 2 ^ 10, maxiter=50 * length(start),derivfree=T),robfunc,varmodel,tau=varmodel$par,...){
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace
	derivfree <- control$derivfree
	
	startlist <- c(start,tau)
	startlist <- as.list(startlist)
	if (is.null(startlist$sigma)){
		dt1 <- c(data,startlist)
    ht <- eval(formula,dt1)
		dv <- as.numeric(ht$predictor) - as.numeric(ht$response)
		#startlist[["sigma"]] <- mad(dv)   # mscale(dv)
	}
	if(derivfree){
		stage1<- optim.NM(objfnc=loss.hetroWM,data=data,start=startlist,formula=formula,varmodel=varmodel,
			robfunc=robfunc,control=control,...)
		method <- fittmethod(methodID= 10,	methodBR=13,subroutine="nl.robhetroWM",lossfunction="loss.hetroWM",subroutineBR="optim.NM")
	}
	else{		
		stage1<- optim.NLM(objfnc=loss.hetroWM,data=data,start=startlist,formula=formula,varmodel=varmodel,
			robfunc=robfunc,control=nlr.control(tolerance=tolerance,maxiter=maxiter),...)
    if(is.Fault(stage1)){
			print("Error In stage 1 Im trying another optimization parameter.")
			stage1<- optim.NM(objfnc=loss.hetroWM,data=data,start=startlist,formula=formula,varmodel=varmodel,
				robfunc=robfunc,control=nlr.control(tolerance=tolerance,maxiter=maxiter),...)
			method <- fittmethod(methodID= 10,	methodBR=13,subroutine="nl.robhetroWM",lossfunction="loss.hetroWM",subroutineBR="optim.NM")
		}
		else{
			method <- fittmethod(methodID= 10,	methodBR=8,subroutine="nl.robhetroWM",lossfunction="loss.hetroWM",subroutineBR="optim.NLM")
		}
	}

	if(is.Fault(stage1)){
		stage1@Fault@FF<-"nl.robhetroWM"
		return(stage1)
	}

	htheta <- stage1$objfnc$value

	## replace latter  vm <- diag(sqrt(as.numeric(stage1$objfnc$varcomp$predictor))/stage2$parameters$sg)
	##rmat <- diag(1.0/sqrt(sqrt(vc)/stage2$parameters$sg))
	vm <- as.numeric(stage1$objfnc$varcomp$predictor)

	if(any(vm<0)) return(Fault(FL="T",FN=19,FF="nl.robhetroWM"))
	rm=diag(1.0/sqrt(vm))
	vm <- diag(vm)
	thetahat <- stage1$parameters[names(formula$par)]
	tauhat <- stage1$parameters[names(varmodel$par)]
	ssq <- sum(stage1$objfnc$ri ^ 2)
  sigma <- ssq / (length(as.numeric(stage1$objfnc$fmod$response)) - formula$p)
	thetahat["sigma"] <- sigma     # delete this later, added scale

	.temphis <- stage1$history
	colnames(.temphis)[2] <- "objfnc"
	result <- nl.fitt.rgn(
      parameters =   thetahat,
      scale =        sigma,
			correlation =  stage1$objfnc$correlation,
			form =         formula,
			response =     stage1$objfnc$fmod$response,
			predictor =    stage1$objfnc$fmod$predictor,
			curvature =    NULL,
			history =      .temphis,
			method =       method,
			data =         as.list(data),
			sourcefnc =    match.call(),
			Fault =        Fault(),
			htheta =       stage1$objfnc$value,
			rho =          stage1$objfnc$rho,
			ri =           stage1$objfnc$ri,
			curvrob =      NULL,
			robform =      robfunc,
			vm =           vm,
			rm=            rm                 # automatically transform will be don by object#
			                                  # gresponse, gpredictor                        #
	)

	vmdata <- as.list(NULL)

	if(is.null( data[[varmodel$independent]])) vmdata[[varmodel$independent]] <- as.numeric( stage1$objfnc$fmod$predictor)
	else vmdata[[varmodel$independent]] <-  data[[varmodel$independent]]
	vmdata[[varmodel$dependent]] <- stage1$objfnc$varcomp$response

  result@hetro<-nl.fitt.rob(
								parameters=    tauhat,
                scale =        sigma,
								form=          varmodel,
								predictor =    stage1$objfnc$varcomp$predictor, 
								response =     stage1$objfnc$ri,
								data =         vmdata,
								method =       NULL,
								others =       list(refvar=stage1$objfnc$refvar)
                            )
	return(result)
}

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.robhetroWM                              |
#|                                                                                 |
#|                              Nov 2012                                           |
#|                                                                                 |
#|                    Hossein Riazoshams, Sweden                                   |
#|                                                                                 |
#+#################################################################################+



