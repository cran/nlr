
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function nlsnmM, optimize function by Nelder Mead.                   | *
#* |                                                                        | *
#* |    argumnts:                                                           | *
#* |      objfnc: any objective function for minimizing, it must contains:  | *
#* |              formula, data and start, extra will be defined in (...)   | *
#* |              the output of objfnc must contains:                       | *
#* |                $value(attr,gradient,hessian), $angmat,$angvec          | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      theta:   starting values, it must contains tau elements           | *
#* |      delta:   for attaining (p+1) initial values.                      | *
#* |         Note: it depends on scale so its better to assigne by program  | *
#* |                                                                        | *
#* |                                                                        | *
#* |      ...:     can be entries objfnc                                    | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
nlsnm<-function(formula,data,start=getInitial(formula,data),delta=NULL,
	control=nlr.control(tolerance=1e-4, maxiter=100 * length(start)),vm=NULL,rm=NULL,...)#eiginv(t(chol(vm))),...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	optfitt <- optim.NM(loss.SSQ,data=data,start=start,formula=formula,vm=vm,rm=rm,control=control,delta=delta,...)
	if(is.Fault(optfitt)) return(optfitt)
	else Fault2 = Fault()
	p <- formula$p
	n<-length(data[[ formula$independent[1] ]])
	dimnames(optfitt$history)[[2]][2]="objfnc"
  sigma <- sqrt(optfitt$history[length(optfitt$history[,1]),"objfnc"]/(n-p))
  names(sigma) = NULL
	optfitt$parameters["sigma"]<- sigma
	ybar <- mean(  as.numeric(optfitt$objfnc$fmod$predictor) )
	ssto <- sum( (as.numeric(optfitt$objfnc$fmod$predictor ) - ybar) ^ 2 )
	if( is.null(vm))
		result <- nl.fitt(
            parameters =  optfitt$parameters,
            scale = sigma,
						correlation =  optfitt$objfnc$correlation,
						form =         formula,
						response =     optfitt$objfnc$fmod$response,
						predictor =    optfitt$objfnc$fmod$predictor, 
#						curvature =    curvatures,
						history =      optfitt$history,
						method =       fittmethod(methodID=			3,
										                  methodBR=			13,
													          	detailBR=			"Nelder Mead",
													          	subroutine=		"nlsnm"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2
					)
	else
		result <- nl.fitt.gn(parameters =  optfitt$parameters,
		        scale = sigma,
						correlation =  optfitt$objfnc$correlation,
						form =         formula,
						response =     transforminv( optfitt$objfnc$fmod$response,rm),
						predictor =    transforminv(optfitt$objfnc$fmod$predictor,rm), 
#						curvature =    curvatures,
						history =      optfitt$history,
						method =    	 fittmethod(methodID=			11,
														methodBR=			13,
														detailBR=			"Nelder Mead",
														subroutine=		"nlsnm"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2,
						vm =           vm,
						rm =           rm,
						gresponse =    optfitt$objfnc$fmod$response,
						gpredictor =   optfitt$objfnc$fmod$predictor
					)

}