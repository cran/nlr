
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nlmest', M-estimate of a nonlinear function.               | *
#* |        Levenberg Marquardt Method.                                     | *
#* |  Note: becarefull to using this function when there is not outlier, it | *
#* |    may not work witout outlier, in this case better to use nlmest      | *
#* |  the problem is in part of two (p2) in hessian its big here.           | *
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

nlmest.LM<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=0.001, minlanda=1 / 2 ^ 10, maxiter=25 * length(start),robscale=T),vm=NULL,rm=eiginv(t(chol(vm))),...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	.tmp <-  (is.null(start$sigma))
	if (.tmp){
	  dt1 <- c(data,start)
	  names(dt1) <- c(names(data),names(start))
	  ht <- eval(formula,dt1)
	  if (is.Fault(ht)) return(ht)
	  dv <- as.numeric(ht$predictor) - as.numeric(ht$response)
	  if(! is.null(vm)) dv <- rm %*% dv		
	  start[["sigma"]] <- mad(dv)   # mscale(dv)
	}
	th <- start											##  with sigma
	theta <- unlist(start[names(start)!="sigma"])	##  without sigma
	theta1 <- theta										##  without sigma
	datalist<-as.list(data)
	.parameters <- names(start)
  p <- length(theta)
	n <- length(data[[1]])
	datalist[.parameters] <- start[.parameters]      #*****datalist contains both parameter vectors and data values
	eol <- F
	iterate <- 0
	iterhist <- NULL
	if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,control=control,...)		## th: with sigma
	else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
	if(any(is.missing(ht$ri))) {
		 return(nl.fitt.rob(Fault=Fault(FN=14,FF="nlmest.LM")))}
	tmp <- as.numeric(ht$ri)
	if(control$robscale) sigma <- mscale(tmp)  #nl.mscale(tmp,robfunc,...)
	else sigma <- sum(abs(tmp)) / n
	if(is.Fault(sigma)) return(sigma)
#lines(data$xr,ht$fmod$predictor)
	if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
	yresp <- as.numeric(ht$fmod$response)
	ybar <- mean(yresp)
	lambda <- 1
						#///*************************************************
	while (!eol)		#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
		grh <- attr(ht$htheta,"gradient")									## g(theta)  grad
		hsh <- attr(ht$htheta,"hessian")									## H(theta)  hess
		ilev <-0
								###################################################################
								###      LM part
								###################################################################
												heig <-eigen(hsh,symmetric=T)			## eig(f+land I) = (eig(F) +lnd)
												lambda <- max(diag(hsh))*2 #heig$values)) * 2	## add the smallest minus eig to
																							## to make all possitive

												repeat{
													ilev<- ilev + 1
													I <- diag(dim(hsh) [2])
													zeta <- hsh + lambda * diag(abs(diag(hsh)))
													zetainv <- eiginv(zeta,stp=F,symmetric=T)
													if(is.Fault(zetainv))  return(nl.fitt.rob(Fault=zetainv))
													delta2 <- zetainv %*% grh
													theta2 <- theta1 - delta2			#### without sigma
													th2 <- as.list(theta2)
													names(th2) <- names(formula$par)
													th2["sigma"]  <- sigma
													if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
													else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
													if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2$Fault))
													cnv <- sum((theta1-theta2)^2)
													if(as.numeric(ht2$htheta) <= as.numeric(ht$htheta) || cnv < tolerance/10^10)
													{
														theta1 <- theta2
														ht <- ht2
														lamnbda <- lambda/10
														break
													}
													else {
														lambda <- lambda * 10
														if( lambda > 1e12) return(nl.fitt.rob(Fault=Fault(FN=14,FF="nlmest.LM")))
													}
												}
		g2 <- attr(ht$rho,"gradient")									## V(heta) = [rho.(r1/sg)   .....  rho.(rn/sg)]T
		tmp <- as.numeric(ht$ri)
		if(control$robscale) sigma <- mscale(tmp)  # nl.mscale(tmp,robfunc,...)
		else sigma <- sum(abs(tmp)) / n
		if(is.Fault(sigma)) return(sigma)
		.expr1 <- t(attr(ht$ri,"gradient")) %*% attr(ht$ri,"gradient")	## J' J  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)
		if(is.Fault(.expr2)) return(nl.fitt.rob(Fault=.expr2)) 		## Note: eiginv only shoud return back Fault
		.expr3 <- attr(ht$ri,"gradient") %*% .expr2						## J (J' J)^-1		n*p
		.expr4 <- .expr3 %*% t(attr(ht$ri,"gradient"))					## J (J' J)^-1 J'  n*n
		.expr5 <- .expr4 %*% g2												## VT = J (J' J)^-1 J' V 	n*1
		angle <- t(g2) %*% .expr5											## V' * VT		1*1
		angle <- angle / sqrt( sum(g2)^2 * sum(.expr5^2) )
		th <-as.list(theta1)													##  without sigma
		names(th) <- names(formula$par)
		th["sigma"]  <- sigma												##  renew with sigma
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
											unlist(th),converge = angle	,ilev=ilev))
		if(angle < tolerance) eol <- T 
		else if(iterate > maxiter) {
			eol <- T
			Fault2 <- Fault(FN=1,FF = "nlmest")
		}
		else {
			if(is.null(vm)) ht <- robloss(formula,data,th,robfunc,control=control,...)					## h(theta) = sum(  rho(ri/sigma)  )
			else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
#			lines(data$xr,ht$fmod$predictor)
			if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
		}
	}		#\\\*********     End Iteration      ****************
			#\\\*************************************************
	if(! is.null(vm)){
		rinv <- eiginv(rm)
		yresp <- rinv %*% as.numeric(ht$fmod$response)
		ypred <- rinv %*% as.numeric(ht$fmod$predictor)
	}
	else{
		 yresp <- as.numeric(ht$fmod$response)
		ypred <- as.numeric(ht$fmod$predictor)
	}
	ybar <- mean(yresp)	
	nlrho <- 1 - sum( ( yresp - ypred )  ^ 2 ) / 
					sum( (ypred - ybar ) ^ 2 )
	curv1 <- curvature(gradient = attr(ht$fmod$predictor,"gradient"),
						hessian = attr(ht$fmod$predictor,"hessian"),
						sigma = th[["sigma"]])
#	curv2 <- curvature(gradient = attr(ht$fmod$predictor,"gradient"),
#						hessian = attr(ht$fmod$predictor,"hessian"),
#						sigma = th[["sigma"]])						
	if(is.null(vm))
		result <- nl.fitt.rob(parameters =  th,
            scale = sigma,
						correlation =  nlrho,
						form =         formula,

						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor,
						curvature =    curv1,
						history =      iterhist,
						method = 		 fittmethod(methodID=			7,
														methodBR=			switch(control$robscale,8,9),       ### rob or nopnrob variance
														detailBR=			"Lev-Marq",
														subroutine=		"nlmest.LM"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =       robfunc
					)
	else
		result <- nl.fitt.rgn(parameters =  th,
		        scale = sigma,
						correlation =  nlrho,
						form =         formula,

						response =     transforminv(ht$fmod$response,rm),
						predictor =    transforminv(ht$fmod$predictor,rm), 
						curvature =    curv1,
						history =      iterhist,
						method = 		 fittmethod(methodID=			7,
														methodBR=			switch(control$robscale,8,9),       ### rob or nopnrob variance
														detailBR=			"Lev-Marq",
														subroutine=		"nlmest.LM"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =      robfunc,
						vm =           vm,
						rm=            rm,
						gresponse =    ht$fmod$response,
						gpredictor =   ht$fmod$predictor
					)
		
	return(result)
}
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of 'nlmest.NLM'                            | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, UPM, INSPEM                        | *
#* |                                                                        | *
#* |                Writen Sep 2008, revised Mar 2009                       | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************


