#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nlmest', M-estimate of a nonlinear function.               | *
#* |       Using combined Gauss-Newton  and Levenberg Marquardt Method.     | *
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
#* |      vm:      variance matrix for generalized minimization if NULL not | *
#* |               generalized.                                             | *
#* |   fixed variance (modified 31/7/2011                                   | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************

nlmest.NLMf<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=0.001, minlanda=1 / 2 ^ 25, maxiter=25 * length(start),robscale=T),vm=NULL,rm=eiginv(t(chol(vm))),...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	if(trace) plot(data$xr,data$yr)
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	.tmp <-  (is.null(start$sigma))
	if(.tmp) return(nl.fitt.rob(Fault=Fault(FL=T,FN=8,FF="nmlest.NLMf")))
	th <- start															##  with sigma
	theta <- unlist(start[names(start)!="sigma"])					##  without sigma
	theta1 <- theta														##  without sigma
	datalist<-as.list(data)
	p <- length(theta)
	n <- length(data[[1]])
	switch(mode(start), 
		"numeric"={.parameters <- names(start)
					start <- as.list(start)
					names(start) <- .parameters; },
		 "list"={.parameters <- names(start); },
		 "NULL"=.parameters <- parameter.names(formula, datalist),
		 return(new("nl.fitt.rob",Fault=Fault(FN=10,FF="nlmest.NLMf")))
	)
	if (!length(.parameters))  return(new("nl.fitt.rob",Fault=Fault(FN=11,FF="nlmest")))
	datalist[.parameters] <- start[.parameters]      #*****datalist contains both parameter vectors and data values
	eol <- F
	iterate <- 0
	iterhist <- NULL
	if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,control=control,...)		## th: with sigma
	else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
#	cat("\n begining ... ht=",as.numeric(ht$htheta),"\n")
	tmp <- as.numeric(ht$ri)
	sigma <- start$sigma
	if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
	if(trace) lines(data$xr,as.numeric(ht$fmod$predictor))
	yresp <- as.numeric(ht$fmod$response)
	ybar <- mean(yresp)
	landa <- 1
	fact <- 2
						#///*************************************************
	while (!eol)		#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
		grh <- attr(ht$htheta,"gradient")									## g(theta)  grad
		hsh <- attr(ht$htheta,"hessian")									## H(theta)  hess
		hshinv <- eiginv(hsh,symmetric=T,stp=F)
		ilev <-0
		if(is.Fault(hshinv))  				#####################################  LM modified
		{
												heig <-eigen(hsh,symmetric=T)			## eig(f+land I) = (eig(F) +lnd)
												lambda <- abs(min(heig$values))*2 		## add the smallest minus eig to
																							## to make all possitive
												repeat{
													ilev<- ilev + 1
													I <- diag(dim(hsh) [2])
													zeta <- hsh + lambda * I#diag((diag(hsh)))
													zetainv <- eiginv(zeta,stp=F,symmetric=T)
#											cat("in singular.....theta 1==",theta1,"\n","values 1=", as.numeric(ht$value))
#											print(is.Fault(zetainv))
													if(is.Fault(zetainv)) return(nl.fitt.rob(Fault=zetainv))
													delta2 <- zetainv %*% grh
													theta2 <- theta1 - delta2			#### without sigma
													th2 <- as.list(theta2)
													th2["sigma"]  <- sigma
													if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
													else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))
#						cat("\n \n \n       in levmarq.....difference... tolerance",diference ,tolerance/10^10,"\n")
#						cat("       .....theta 2...theta1==",theta1,theta2,"\n"," \n values 2..1=",as.numeric(ht2$value), as.numeric(ht$value))
#						cat("\n     .....lambda... delta",lambda,delta2 ," \n")
													cnv <- sum((theta1-theta2)^2)
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))
													if( (as.numeric(ht2$htheta) <= as.numeric(ht$htheta)) || (cnv < tolerance/10^10) )
													{
														theta1 <- theta2
														ht <- ht2
														lambda<-lambda/10
														break

													}
													else {
														lambda <- lambda * 10
														if( lambda > 1e12) return(nl.fitt.rob(Fault=Fault(FN=14,FF="nlmest.NLM")))
													}
												}  #####################################  enf of LM
		}
							#####################     end of singular case                 ####################
							###################################################################################
		else {  			#####################   Possitive definit case  ###################################
						delta1 <- hshinv %*% grh
						theta2 <- theta1 - landa*delta1						
						th2 <- as.list(theta2)
						th2["sigma"] <- sigma
						if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
						else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
						cnvg <- F
						diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))
#				cat("positive definit case","iterate=",iterate, "difference tol=",diference ,tolerance/1000,"\n")
#				cat("theta 2...theta1==",theta1,theta2,"\n","values 2..1=",as.numeric(ht2$htheta), as.numeric(ht$htheta))
#				cat("\n landa... delta",landa,delta1,landa*delta1,"\n  \n")
						temp1 <- as.numeric(ht2$value)
						temp2 <-  as.numeric(ht$value)
						if(as.numeric(ht2$htheta) > as.numeric(ht$htheta)){
																	###############################
							while(landa >= minlanda){ 			####### iteration landa     ###
									landa <- landa / fact
									theta2 <- theta1 - landa * delta1
									th2 <- as.list(theta2)
									th2["sigma"] <- sigma
									if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
									else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
									diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))
#						cat("\n \n \n       in landa.....difference... tolerance",diference ,tolerance/1000,"\n")
#						cat("       .....theta 2...theta1==",theta1,theta2,"\n","values 2..1=",as.numeric(ht2$value), as.numeric(ht$value))
#						cat("\n     .....landa... delta",landa,delta1,landa*delta1,"\n \n")
									if(as.numeric(ht2$htheta) <= as.numeric(ht$htheta))
										{
											theta1 <- theta2
											ht <- ht2
											cnvg <- T
											break
										}
							}										#######  iteration landa    ###
																	###############################
							landa <- min(1,fact*landa)
							#landa <- 1
						}
						else{
							theta1 <- theta2
							ht <- ht2
							#landa<-landa*fact
							cnvg <- T
						}
						if(! cnvg) {						# ***********************************************
															# ****  non convergance using LM again.      ****
														#	 return(nlmest.LM(formula,data=data,start=th2,robfunc=robfunc,control=control,rm=rm,...))
												heig <-eigen(hsh,symmetric=T)			## eig(f+land I) = (eig(F) +lnd)
												lambda <- abs(min(heig$values)) * 2	## add the smallest minus eig to
																							## to make all possitive
												repeat{
													ilev<- ilev + 1
													I <- diag(dim(hsh) [2]) 			 ### this  ones better faster converge.
													zeta <- hsh + lambda * I
													zetainv <- eiginv(zeta,stp=F,symmetric=T)
													if(is.Fault(zetainv)) return(nl.fitt.rob(Fault=zetainv))
													delta2 <- zetainv %*% grh
													theta2 <- theta1 - delta2			#### without sigma
													th2 <- as.list(theta2)
													th2["sigma"]  <- sigma
													if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
													else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
													cnv <- sum((theta1-theta2)^2)
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))
													if(as.numeric(ht2$htheta) <= as.numeric(ht$htheta))												
													{
														theta1 <- theta2
														ht <- ht2
														break
													}
													else {
														lambda <- lambda * 10
														#if( lambda > 1e12) return(nl.fitt.rob(Fault=Fault(FN=14,FF="nlmest.NLM")))
														if( lambda > 1e12) return(nlmest.LM(formula,data=data,start=th2,robfunc=robfunc,control=control,vm=vm,rm=rm,...))
													}
												}
						}									# *******      enad of LM again             ********
						  									# **************************************************
		}			#######################     end positive                                      ############
					#######################  from above theta1 & ht must be returned back         ############
					##########################################################################################
		g2 <- attr(ht$rho,"gradient")									## V(heta) = [rho.(r1/sg)   .....  rho.(rn/sg)]T
		tmp <- as.numeric(ht$ri)
		.expr1 <- t(attr(ht$ri,"gradient")) %*% attr(ht$ri,"gradient")	## J' J  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)
		if(is.Fault(.expr2)){
			heig <-eigen(.expr1,symmetric=T)
			lambda <- abs(min(heig$values)) * 2
			I <- diag(dim(.expr1) [2])
			.expr1 <- .expr1 + lambda * I
			.expr2 <- eiginv(.expr1,symmetric=T,stp=F)			
		}
		if(is.Fault(.expr2)) return(nl.fitt.rob(Fault=.expr2))
		.expr3 <- attr(ht$ri,"gradient") %*% .expr2						## J (J' J)^-1		n*p
		.expr4 <- .expr3 %*% t(attr(ht$ri,"gradient"))					## J (J' J)^-1 J'  n*n
		.expr5 <- .expr4 %*% g2												## VT = J (J' J)^-1 J' V 	n*1
		angle <- t(g2) %*% .expr5											## V' * VT		1*1
		angle <- angle / sqrt( sum(g2^2) * sum(.expr5^2) )
		th <-as.list(theta1)													##  without sigma
		th["sigma"]  <- sigma												##  renew with sigma
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
											unlist(th),converge = angle	,ilev=ilev))
		if(angle < tolerance || diference < tolerance/100) eol <- T 
		else if(iterate > maxiter) {
			eol <- T
			Fault2 <- Fault(FN=1,FF = "nlmest")
		}
		else {
			if(is.null(vm)) ht <- robloss(formula,data,th,robfunc,control=control,...)  ## h(theta) = sum(  rho(ri/sigma)  )
			else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
			if(trace) lines(data$xr,ht$fmod$predictor)
			if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
			tmp <- as.numeric(ht$ri)
			th["sigma"]  <- sigma												##  renew with sigma				
										}
			#\\\*********                        ****************
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
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor,
						curvature =    curv1,
						history =      iterhist,
						method = 				fittmethod(methodID=			7,
															methodBR=			3,       ### fixed variance
															detailBR=			"Newton Lev-Marq",
															subroutine=		"nlmest.NLMf"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =      robfunc
					)
	else

		result <- nl.fitt.rgn(parameters =  th,
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor, 
						curvature =    curv1,
						history =      iterhist,
						method = 			fittmethod(methodID=			8,
															methodBR=			3,       ### fixed variance
															detailBR=			"Newton Lev-Marq",
															subroutine=		"nlmest.NLMf"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =       robfunc,
						vm =           vm,
						rm=            rm,
						gresponse =    transform(ht$fmod$response,rm),
						gpredictor =   transform(ht$fmod$predictor,rm)
					)
	return(result)
}
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of optim.NLM'                              | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, UPM, INSPEM                        | *
#* |                                                                        | *
#* |                Recoded from 'nlmest' in Jan 2010                       | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************









