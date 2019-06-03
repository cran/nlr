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
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************

nlmest.NLM<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=0.0001, minlanda=1 / 2 ^ 25, maxiter=25 * length(start)),vm=NULL,rm=NULL,...)
{
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace
  if(is.null(rm) && ! is.null(vm)) rm <- eiginv(t(chol(vm)))
	#if(is.null(vm)) rm <- NULL
	#else rm <- eiginv(t(chol(vm)))
	loc.start <- start
	if(trace) plot(data[[formula$independent]],data[[formula$dependent]])
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	if (is.null(loc.start$sigma)){
		dt1 <- c(data,loc.start)
    names(dt1) <- c(names(data),names(loc.start))
		ht <- eval(formula,dt1)
    if (is.Fault(ht)) return(ht)
		dv <- as.numeric(ht$predictor) - as.numeric(ht$response)
		if(! is.null(vm)) dv <- rm %*% dv		
		loc.start[["sigma"]] <- mad(dv)   # mscale(dv)
	}
	th <- loc.start                                               ##  with sigma
	theta <- unlist(loc.start[names(loc.start)!="sigma"])         ##  without sigma
	theta1 <- theta                                               ##  without sigma
	datalist<-as.list(data)
  .parameters <- names(start)
	p <- length(theta)
	n <- length(data[[1]])
	datalist[.parameters] <- loc.start[.parameters]      #*****datalist contains both parameter vectors and data values
	eol <- F
	iterate <- 0
	iterhist <- NULL
	lambda <-1
	if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,control=control,...)		## th: with sigma
	else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
  if(is.Fault(ht)) return(nl.fitt.rob(Fault=ht))		
	tmp <- as.numeric(ht$ri)
	if(control$robscale) sigma<-nl.mscale(tmp,robfunc,...) # <- mad(tmp)
	else sigma  <- sum(abs(tmp)) / n
	if(is.Fault(sigma)) return(sigma)

	if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
	if(trace) lines(data[[formula$independent]],as.numeric(ht$fmod$predictor))

	yresp <- as.numeric(ht$fmod$response)
	ybar <- mean(yresp)
	landa <- 1
	alpha <- 1
	fact <- 2
						#///*************************************************
	while (!eol)		#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
#   cat("iterate ============================================= ",iterate,"\n")
		grh <- attr(ht$htheta,"gradient")									## g(theta)  grad
		hsh <- attr(ht$htheta,"hessian")									## H(theta)  hess
    hshinv <- eiginv(hsh,symmetric=T,stp=F)
		ilev <-0
		if(is.Fault(hshinv))  				#####################################  LM modified
		{
												heig <-eigen(hsh,symmetric=T)			## eig(f+land I) = (eig(F) +lnd)
												#lambda <-  abs(min(heig$values))*2 		## add the smallest minus eig to
												lambda <-1
												damping <-  diag(dim(hsh) [2])
												.temp1 <-heig$values
												.temp1[.temp1<0]<-(1e-8)-.temp1[.temp1<0]
												.temp1 <- diag(.temp1)
												.temp2<-(heig$vectors) %*% .temp1
												zeta <- .temp2 %*% t(heig$vectors) 
												zetainv <- indifinv(zeta,stp=F,symmetric=T)
#                       print(zeta)
												if(is.Fault(zetainv)) return(nl.fitt.rob(Fault=zetainv))
												delta2 <- zetainv %*% grh
												if(is.matrix(delta2)) delta2<-as.numeric(delta2)
												repeat{
													ilev<- ilev + 1
													theta2 <- theta1 - lambda * delta2			#### without sigma
#                         cat("gradient =",grh>0,"\n")
#                         cat("lambda=",lambda,"delta2=",delta2,"increament = ",- lambda * delta2	,"\n")
#                         cat("theta 2 = ",theta2,"theta1 = ", theta1,"\n")
													th2 <- as.list(theta2)
													names(th2) <- names(formula$par)
													th2["sigma"]  <- sigma
													if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
													else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
													if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
													cnv <- sum((theta1-theta2)^2)/(sum(theta1^2)+(1e-8))
                          if( (as.numeric(ht2$htheta) <= as.numeric(ht$htheta)))# || (cnv < tolerance/10000) )
													{
														theta1 <- theta2
														ht <- ht2
														#alpha <- alpha * fact
														#if(alpha > 1e10) alpha <- 1
														break
													}
													else {
														lambda <- lambda / fact
														if( lambda <minlanda){
															#sigma<- mscale(ht2$ri) 
															#lambda <- lambda * 20 
															return(nl.fitt.rob(Fault=Fault(FN=14,FF="nlmest.NLM")))
														}
														#alpha <- alpha / fact
														#if(alpha < minlanda) alpha <- 1
													}
												}  #####################################  enf of LM
		}
							#####################     end of singular case                 ####################
							###################################################################################
		else {  			#####################   Possitive definit case  ###################################
						landa <- 1
						delta1 <- hshinv %*% grh
						theta2 <- theta1 - landa*delta1						
						th2 <- as.list(theta2)
						names(th2) <- names(formula$par)
						th2["sigma"] <- sigma
#           cat("\n positive definit case---------------------------------------------------- \n")
						if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
						else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
						if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
						cnvg <- F
						diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
						gradpk <- t(grh) %*% delta1
						rate <- 1e-4 * landa * gradpk
#						cat("\n iterate=",iterate, "difference tol=",diference ,tolerance/1000,"\n")
#						cat("values 2..1=",as.numeric(ht2$htheta),as.numeric(ht$htheta))
#						cat("\n theta 1= ",theta1,"theta2 = ",theta2)
#						cat("\n landa... delta",landa,delta1,"theta 2=",theta2)
#						cat("\n sigma=",sigma,"\n")
						.temp1<-as.numeric(ht2$htheta)
						.temp2<-as.numeric(ht$htheta)
						if(	is.nan(.temp1) || is.infinite(.temp1) ||
								is.nan(.temp1) || is.infinite(.temp1)) return(nl.fitt.rob(Fault=Fault(FN=18,FF="nlmest.NLM")))
						if(as.numeric(ht2$htheta) > as.numeric(ht$htheta)){
																###############################
							while(landa >= minlanda){ 			####### iteration landa     ###
									theta2 <- theta1 - landa * delta1
									th2 <- as.list(theta2)
									names(th2) <- names(formula$par)
									th2["sigma"] <- sigma
									if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
									else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
									if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
									diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
									if(as.numeric(ht2$htheta) <= as.numeric(ht$htheta))
										{
											theta1 <- theta2
											ht <- ht2
											cnvg <- T
											break
										}
									landa <- landa / fact
							}										#######  iteration landa    ###
																	###############################
							landa <- min(1,fact*landa)
							#landa <- 1
						}
						else{
							theta1 <- theta2
							ht <- ht2
							landa<-landa*fact
							cnvg <- T
						}#######  iteration landa    ###
						###############################
						if(!cnvg){
															# ***********************************************
															# ****  non convergance using LM again.      ****
															# ****                                       ****
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
													names(th2) <- names(formula$par)
													th2["sigma"]  <- sigma
													if(is.null(vm)) ht2 <- robloss(formula,data,th2,robfunc,control=control,...)
													else ht2 <- robloss.gn(formula,data,th2,robfunc,rm,control=control,...)
													if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
														
													cnv <- sum((theta1-theta2)^2)
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/abs(as.numeric(ht$htheta))
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
		if (control$robscale) sigma <-nl.mscale(tmp,robfunc,...)# mscale(tmp)
		else sigma <- sum(abs(tmp)) / n
		if(is.Fault(sigma)) return(sigma)			
#		cat("\n sigma = ",sigma,"\n")
		.expr1 <- t(attr(ht$ri,"gradient")) %*% attr(ht$ri,"gradient")	## J' J  p*p
		.expr2 <- indifinv(.expr1,symmetric=T,stp=F)
		if(is.Fault(.expr2)){
		  if(any(is.na(.expr1)) || any(is.inf(.expr1))) return(nl.fitt.rob(Fault = Fault(FN=4,FF="nlmest.NLM")))
			heig <-eigen(.expr1,symmetric=T)
			lambda <- abs(min(heig$values)) + 1.01
			I <- diag(dim(.expr1) [2])
			.expr1 <- .expr1 + lambda * I
			.expr2 <- indifinv(.expr1,symmetric=T,stp=F)			
			if(is.Fault(.expr2)) .expr2 <-ginv(.expr1)
		}
		if(is.Fault(.expr2)) return(nl.fitt.rob(Fault=.expr2))
		.expr3 <- attr(ht$ri,"gradient") %*% .expr2						## J (J' J)^-1		n*p
		.expr4 <- .expr3 %*% t(attr(ht$ri,"gradient"))					## J (J' J)^-1 J'  n*n
		.expr5 <- .expr4 %*% g2												## VT = J (J' J)^-1 J' V 	n*1
		angle <- t(g2) %*% .expr5											## V' * VT		1*1
		angle <- angle / sqrt( sum(g2^2) * sum(.expr5^2) ) 
    angle <- abs(angle)
		th <-as.list(theta1)												##  without sigma
		names(th) <- names(formula$par)
		th["sigma"]  <- sigma												##  renew with sigma
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
											unlist(th),converge = angle	,ilev=ilev))
#   cat("1111 angle 2= ", angle,"diference =",diference ,"\n")
		if(is.missing(angle)||	is.nan(angle) || is.inf(angle)) return(nl.fitt.rob(Fault=Fault(FN=16,FF="nlmest.NLM")))
    if(angle < tolerance || diference < tolerance/1000) eol <- T 
		else if(iterate > maxiter) {
			eol <- T
			Fault2 <- Fault(FN=1,FF = "nlmest")
		}
		else {
			if(is.null(vm)) ht <- robloss(formula,data,th,robfunc,control=control,...)  ## h(theta) = sum(  rho(ri/sigma)  )
			else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
			if(trace) lines(data[[formula$independent]],ht$fmod$predictor)
#			cat("222222 sigma= ",sigma,"\n")
			if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
			tmp <- as.numeric(ht$ri)
			if (control$robscale) sigma <-nl.mscale(tmp,robfunc,...)   # mscale(tmp)
			else sigma <- sum(abs(tmp)) / n
			th["sigma"]  <- sigma												##  renew with sigma				
										}
			#\\\*********                        **********************************************************************************
	}		#\\\*********     End Iteration      **********************************************************************************
			#\\\*******************************************************************************************************************
	#//////////*******************************************************************************************************************
	#//////////*******************************************************************************************************************
	#//////////********     preparing output   ***********************************************************************************
	#//////////********                        ***********************************************************************************
	if(! is.null(vm)){
		rinv <- indifinv(rm)
	  #rinv <- solve(rm)
		#yresp <- rinv %*% as.numeric(ht$fmod$response)
		#ypred <- rinv %*% as.numeric(ht$fmod$predictor)
		yresp <- as.numeric(ht$fmod$response)
		ypred <- as.numeric(ht$fmod$predictor)		
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
mtbr <- switch(as.numeric(control$robscale)+1,5,4)

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
														methodBR=			switch(as.numeric(control$robscale)+1,5,4),       ### rob or nopnrob variance
														detailBR=			"Newton Lev-Marq",
														subroutine=		"nlmest.NLM"),
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
            scale = sigma,
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor, 
						curvature =    curv1,
						history =      iterhist,
						method = 	  	 fittmethod(methodID=			8,
									              		methodBR=			switch(as.numeric(control$robscale)+1,5,4),       ### rob or nopnrob variance
							              				detailBR=			"Newton Lev-Marq",
									              		subroutine=			"nlmest.NLM"),
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
						gresponse =    transformNR(ht$fmod$response,rm),
						gpredictor =   transformNR(ht$fmod$predictor,rm)
					)
  if(is.null(vm)) 	result@method@lossfunction <-"robloss"
  else result@method@lossfunction <- "robloss.gn"
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









