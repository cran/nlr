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

dfrmest.NLM<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=0.01, minlanda=1 / 2 ^ 25, maxiter=25 * length(start)),vm=NULL,rm=NULL,...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	if(is.null(vm)) rm <- NULL
	else rm <- eiginv(t(chol(vm)))
	loc.start <- start
	if(trace) plot(data[[formula$independent]],data[[formula$dependent]])
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	if (is.null(loc.start$sigma)){
		dt1 <- c(data,loc.start)
		ht <- eval(formula,dt1)
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
  names(th) <- names(loc.start) 
  ht <- dfr.robloss(formula,data,th,robfunc,control=control,rmat=rm,...)		## th: with sigma
	if(is.Fault(ht)) return(nl.fitt.rob(Fault=ht))
		
	tmp <- as.numeric(ht$ri)
	if(control$robscale) sigma<-mscale(tmp) # <- nl.mscale(tmp,robfunc,...)
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
		#cat("iterate ============================================= ",iterate,"\n")
		grh <- attr(ht$htheta,"gradient")									## g(theta)  grad
		hsh <- attr(ht$htheta,"hessian")									## H(theta)  hess
		hshinv <- eiginv(hsh,symmetric=T,stp=F)
		ilev <-0
		if(is.Fault(hshinv))  				#####################################  LM modified
		{
#												cat("singularity.. iteration= ",iterate,"\n")
												heig <-eigen(hsh,symmetric=T)			## eig(f+land I) = (eig(F) +lnd)
												#lambda <-  abs(min(heig$values))*2 		## add the smallest minus eig to
												lambda <-1
												damping <-  diag(dim(hsh) [2])
												.temp1 <-heig$values
												.temp1[.temp1<0]<-(1e-8)-.temp1[.temp1<0]
												.temp1 <- diag(.temp1)
												.temp2<-(heig$vectors) %*% .temp1
												zeta <- .temp2 %*% t(heig$vectors) 
												zetainv <- eiginv(zeta,stp=F,symmetric=T)
												if(is.Fault(zetainv)) return(nl.fitt.rob(Fault=zetainv))
												delta2 <- zetainv %*% grh
												repeat{
													ilev<- ilev + 1
													theta2 <- theta1 - lambda * delta2			#### without sigma
													th2 <- as.list(theta2)
													names(th2) <- names(formula$par)
													th2["sigma"]  <- sigma
                          ht2 <- dfr.robloss(formula,data,th2,robfunc,control=control,rmat=rm,...)
#             						cat("as.numeric(ht2$htheta)=",as.numeric(ht2$htheta),"as.numeric(ht$htheta)=",as.numeric(ht$htheta))
													if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
													cnv <- sum((theta1-theta2)^2)/(sum(theta1^2)+(1e-8))
													#cat("lambda = ",lambda,"(as.numeric(ht2$htheta)= ",as.numeric(ht2$htheta) ,as.numeric(ht$htheta),"\n")
													if( (as.numeric(ht2$htheta) <= as.numeric(ht$htheta)))# || (cnv < tolerance/10000) )
													{
														theta1 <- theta2
														ht <- ht2
														#alpha <- alpha * fact
														if(alpha > 1e10) alpha <- 1
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
														if(alpha < minlanda) alpha <- 1
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
#						cat("\n positive definit case---------------------------------------------------- \n")
#						cat("as.numeric(ht$htheta)=",as.numeric(ht$htheta))
#						cat("theta 1"=,theta1,"\n","theta 2 and delta=",theta2,"\n",delta1,"\n")
						ht2 <- dfr.robloss(formula,data,th2,robfunc,control=control,rmat=rm,...)
#						cat("as.numeric(ht2$htheta)=",as.numeric(ht2$htheta),"as.numeric(ht$htheta)=",as.numeric(ht$htheta))
						if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
						cnvg <- F
						diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
						gradpk <- t(grh) %*% delta1
						rate <- 1e-4 * landa * gradpk
						.temp1<-as.numeric(ht2$htheta)
						.temp2<-as.numeric(ht$htheta)
						if(is.missing(.temp1) ||	is.nan(.temp1) || is.inf(.temp1) ||
							is.missing(.temp1) ||	is.nan(.temp1) || is.inf(.temp1)) return(nl.fitt.rob(Fault=Fault(FN=18,FF="nlmest.NLM")))
						if(as.numeric(ht2$htheta) > as.numeric(ht$htheta)){
																###############################
							while(landa >= minlanda){ 			####### iteration landa     ###
									theta2 <- theta1 - landa * delta1
									th2 <- as.list(theta2)
									names(th2) <- names(formula$par)
									th2["sigma"] <- sigma
                  ht2 <- dfr.robloss(formula,data,th2,robfunc,control=control,rmat=rm,...)
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
													ht2 <- dfr.robloss(formula,data,th2,robfunc,control=control,rmat=rm,...)
                          if(is.Fault(ht2)) return(nl.fitt.rob(Fault=ht2))
														
													cnv <- sum((theta1-theta2)^2)
													diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/as.numeric(ht$htheta)
													if(as.numeric(ht2$htheta) <= as.numeric(ht$htheta))												
													{
														theta1 <- theta2
														ht <- ht2
														break
													}
													else {
														lambda <- lambda * 10
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

		if (control$robscale) sigma <-mscale(tmp)   # nl.mscale(tmp,robfunc,...)
		else sigma <- sum(abs(tmp)) / n
		if(is.Fault(sigma)) return(sigma)			
		th <-as.list(theta1)												##  without sigma
		names(th) <- names(formula$par)
    th["sigma"]  <- sigma												##  renew with sigma
		
		.expr1 <- t(attr(ht$ri,"gradient")) %*% attr(ht$ri,"gradient")	## J' J  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)

		if(is.Fault(.expr2)){
			cvg <- (diference < tolerance/100)
			iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
						unlist(th),converge = diference, ilev=ilev))
		}	
		else{
			.expr3 <- attr(ht$ri,"gradient") %*% .expr2						## J (J' J)^-1		n*p
			.expr4 <- .expr3 %*% t(attr(ht$ri,"gradient"))					## J (J' J)^-1 J'  n*n
			.expr5 <- .expr4 %*% g2												## VT = J (J' J)^-1 J' V 	n*1
			angle <- t(g2) %*% .expr5											## V' * VT		1*1
			angle <- angle / sqrt( sum(g2^2) * sum(.expr5^2) ) 
			cvg <- (angle < tolerance || diference < tolerance)
			iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
								unlist(th),converge = angle	,ilev=ilev))
		}
		if(cvg) eol <- T 
		else if(iterate > maxiter) {
			eol <- T
			Fault2 <- Fault(FN=1,FF = "nlmest")
		}
		else {

			ht <- dfr.robloss(formula,data,th,robfunc,control=control,rmat=rm,...)  ## h(theta) = sum(  rho(ri/sigma)  )
			if(trace) lines(data[[formula$independent]],ht$fmod$predictor)
			if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
			tmp <- as.numeric(ht$ri)
			if (control$robscale) sigma <-mscale(tmp)   # nl.mscale(tmp,robfunc,...)
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
		rinv <- eiginv(rm)
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
	if(is.null(vm))
		result <- nl.fitt.rob(parameters =  th,
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor,
						curvature =    NULL,
						history =      iterhist,
						method = 		 fittmethod(methodID=			7,
														methodBR=			switch(control$robscale,4,5),       ### rob or nopnrob variance
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
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor, 
						curvature =    NULL,
						history =      iterhist,
						method = 		 fittmethod(methodID=			8,
														methodBR=			switch(control$robscale,4,5),       ### rob or nopnrob variance
														detailBR=			"Newton Lev-Marq",
														subroutine=		"nlmest.NLM"),
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








