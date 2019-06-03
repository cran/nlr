#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nlmmest', MM-estimate of a nonlinear function.             | *
#* |       Using reweighted  Method, by stromberg.                         | *
#* |    argumnts:                                                           | *
#* |      formula: 'nl.form' object, the function mode.                     | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      start:   starting values, it must contains 'sigma', selstart      | *
#* |               for nl.form object is not created yet, take cre of it.   | *
#* |      robfunc: obust function, it must be functionformat untill know,   | *
#* |               not a nl.form, in feature it should be modfied.          | *
#* |      ...:     can be entries for robust loss function parameters.      | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
nlmest.RWT <- function(formula , data , start = getInitial(formula , data) , robfunc ,
	control=nlr.control(tolerance = 0.001,minlanda = 1 / 2 ^ 25 , maxiter = 25 * length(start),trace=F),vm=NULL,rm=eiginv(t(chol(vm))) ,...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	Fault2 <- Fault()                                                          ###### no error, automatically will be created by the Fault object
	loc.start <- start	
	if (is.null(loc.start$sigma)){
		dt1 <- c(data,loc.start)
    names(dt1) <- c(names(data),names(loc.start))
		ht <- eval(formula,dt1)
		dv <- as.numeric(ht$predictor) - data[[formula$dependent]]
		if(! is.null(vm)) dv <- rm %*% dv
		loc.start[["sigma"]] <- mscale(dv)
	}
	th <- loc.start
	names(th) <- names(formula$par)  
	theta <- unlist(start[names(th)!="sigma"])                              ## without sigma, numeric
	theta1 <- theta                                                            ## without sigma
	datalist <- as.list(data)
	p <- length(formula$p)
	n <- length(data[[2]])
	switch(mode(start) , "numeric" = {
				.parameters <- names(start)
				start <- as.list(start)
				names(start) <- .parameters
			} , "list" = {
				.parameters <- names(start)
			} , "NULL" = .parameters <- parameter.names(formula , datalist) ,
			return(new("nl.fitt.rob",Fault=Fault(FN=10,FF="nlmest.RWT")))
	)
	if (!length(.parameters))  return(new("nl.fitt.rob",Fault=Fault(FN=11,FF="nlmest.RWT")))
	datalist[.parameters] <- loc.start[.parameters]                           ##*****datalist contains both parameter vectors and data values
	eol <- F
	iterate <- 0
	iterhist <- NULL
	if(trace) plot(data[[formula$independent]],data[[formula$dependent]])
  if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,robscale=T,...)		## th: with sigma
	else ht <- robloss.gn(formula,data,th,robfunc,rm,robscale=T,...)
	if(trace) lines(data[[formula$independent]],as.numeric(ht$fmod$predictor))
	z <- as.numeric(ht$ri)
	if(is.Fault(ht)) return(ht)
	ybar <- mean(ht$fmod$response)
	scale <- mscale(z)   ## <- nl.mscale(tmp,robfunc,...)
	scale2 <- scale * robfunc$arguments$k1   
	rb <- ht$rho
	rho <- as.numeric(ht$htheta)
	yresp <- as.numeric(ht$fmod$response)
	#///*************************************************
	while (!eol)#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
		G <- attr(ht$fmod$predictor , "gradient") / robfunc$arguments$k1 / scale
		w <- attr(rb , "weight")
		wG <- sqrt(w) %*% G
		ps <- attr(rb , "gradient")
		d1 <- crossprod(G , ps)
		d2 <- crossprod(wG)
		max2 <- max(d2)
		if(iterate == 1)
		{
			if(all(w == 0)) return(Fault(FN=21))
			lambda <- sqrt(mean(d2 ^ 2)) / 1000
			Imat <- diag(dim(G)[2])
		}
		bold <- theta1
		old <- rho
		ilev <- 0
		repeat
		{
			ilev <- ilev+1
			hs <- d2 + lambda * Imat
			aa <- eiginv(hs,symmetric=T)
			db <- aa %*% d1
			theta1 <- bold + db
			th <- as.list(theta1)
			names(th) <- names(formula$par)
			th["sigma"] <- scale
			if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,robscale=T,...)		## th: with sigma
			else ht <- robloss.gn(formula,data,th,robfunc,rm,robscale=T,...)
			rb <- ht$rho
			rhf <- as.numeric(rb)
			rho <- as.numeric(ht$htheta)
			if(rho < old -1e-15 || rho == 0) break
			if(lambda / max2 > 1e15)
			{
				theta1 <- bold
				return(Fault(FN=14,FF="nlmest.RWT"))
#				warning("Levenberg tolerance not acheived")
				break
			}
			lambda <- 2 * lambda
		}
		if(trace)  plot(data[[formula$independent]],as.numeric(ht$fmod$predictor))
		if(lambda / max2 > 1e15) break
		if(ilev == 1) lambda <- lambda / 10
		converge <- crossprod(d1 , db)		                                       ###test for convergence
		if(converge < 1e-2) break
		if(iterate > maxiter)
		{
			warning("mmnl: max iteration exceed")
			break
		}
		iterhist <- rbind( iterhist , c( iteration = iterate , objfnc = as.numeric(ht$htheta) , unlist(th) , converge = converge,ilev=ilev ) )
		if(converge < tolerance) eol <- T
		else if(iterate > maxiter)
		{ 
			eol <- T
			Fault2 <- Fault(FN = 1 , FF = "nlmest.RWT")
		}
	}
	nlrho <- 1 - sum( ( yresp - as.numeric(ht$fmod$predictor) ) ^ 2 ) / sum( (as.numeric(ht$fmod$predictor) - ybar ) ^ 2 )
	curv1 <- curvature(gradient = attr(ht$fmod$predictor,"gradient"),
						hessian = attr(ht$fmod$predictor,"hessian"),
						sigma = th[["sigma"]])
#	curv2 <- curvature(gradient = attr(ht$fmod$predictor,"gradient"),
#						hessian = attr(ht$fmod$predictor,"hessian"),
#						sigma = th[["sigma"]])						
	if(is.null(vm))
		result <- nl.fitt.rob(
						parameters = 			th , 
						correlation = 		nlrho , 
						form = 				formula , 
						response = 			ht$fmod$response , 
						predictor = 			ht$fmod$predictor , 
						curvature = 			curv1 , 
						history = 				iterhist , 
						method = 				fittmethod(methodID=			7,
															methodBR=			2,
															detailBR=			"Reweighting",
															subroutine=		"nlmest.SD"),
						data = 				as.list(data) , 
						sourcefnc =           match.call(),
						Fault = 				Fault2 , 
						htheta = 				ht$htheta , 
						rho = 					ht$rho , 
						ri = 					ht$ri , 
						curvrob =           NULL,
						robform =				robfunc)
	else
		result <- nl.fitt.rgn(
						parameters = 			th , 
						correlation = 		nlrho , 
						form = 				formula , 
						response = 			ht$fmod$response , 
						predictor = 			ht$fmod$predictor , 
						curvature = 			curv1 , 
						history = 			iterhist , 
						method = 				fittmethod(methodID=			8,
															methodBR=			2,
															detailBR=			"Reweighting",
															subroutine=		"nlmest.SD"),
						data = 				as.list(data) , 
						sourcefnc =           match.call(),
						Fault = 				Fault2 , 
						htheta = 				ht$htheta , 
						rho = 					ht$rho , 
						ri = 					ht$ri , 
						curvrob =           NULL,
						robform =				robfunc,
						vm =                vm,
						rm=                 rm,
						gresponse =         transform(ht$fmod$response,rm),
						gpredictor =        transform(ht$fmod$predictor,rm))
	return(result)
}
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of 'nl_mmest'                              | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, UPM, INSPEM                        | *
#* |                                                                        | *
#* |                Writen Sep 2008, revised OCT 2008                       | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
