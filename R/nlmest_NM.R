#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nlmest.NM', M-estimate of a nonlinear function.            | *
#* |       Using Nelder Mead, derivative free.                              | *
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

nlmest.NM<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=1e-4, minlanda=1 / 2 ^ 25, maxiter= 100 * length(start),robscale=T)
	,vm=NULL,rm=eiginv(t(chol(vm))),delta=NULL,...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace
	datalist<-as.list(data)
	p <- formula$p
	n <- length(data[[1]])
	loc.start <- start
	if(trace) plot(data[[formula$independent]],data[[formula$dependent]])
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	if (is.null(loc.start$sigma)){
		dt1 <- c(data,loc.start)
		ht0 <- eval(formula,dt1)
		dv <- as.numeric(ht0$predictor) - as.numeric(ht0$response)
		if(! is.null(vm)) dv <- rm %*% dv		
		loc.start[["sigma"]] <-  mad(dv) #mscale(dv)  
	}
	th <- loc.start                                               ##  with sigma

	sigma <- loc.start[["sigma"]] 
	b<- unlist(loc.start[names(loc.start)!="sigma"])              ## without sigma
	if(is.null(delta)){
    delta<-b * .1
    delta[b==0] <- 0.1
	}
	theta <- rep(b,p)                                             ## without sigma
	dim(theta) <- c(p,p)
	theta <- t(theta)
	diag(theta) <- diag(theta)+delta
	theta <- rbind(b-delta,theta)                                      ## p+1 rows of parameters, exclude sigma

	th <- theta                                                  ##  include sigma
	colnames(theta) <- names(start)[names(start)!="sigma"]
	th <- cbind(th, rep(loc.start[["sigma"]],p+1))
	colnames(th) <- names(loc.start)

	y <- rep(0,p+1)
	iterhist <- NULL
	ht <- list(NULL)

  for (i in 1:(p+1)){
		dataloss <- th[i,]

		if(is.null(vm)) 	ht[[i]] <- robloss(formula,data,start=dataloss ,robfunc,control=control,...)		## th: with sigma
		else ht[[i]] <- robloss.gn(formula,data,dataloss ,robfunc,rm,control=control,...)

		if(is.Faultwarn(ht[[i]])){
			ht[[i]]@FN <- "nlmest.NL"
			return(ht[[i]])
		}
		y[i] <- as.numeric(ht[[i]]$htheta)
		iterhist <- rbind(iterhist,c(iteration = 0,objfnc = as.numeric(ht[[i]]$htheta),dataloss,converge = 0))	
		if(trace) lines(data[[formula$independent]],as.numeric(ht[[i]]$fmod$predictor))
	}
	eol <- F
	iterate <- 0
	psum <- apply(theta,2,sum)                                        ###  sum or mean   ????? corect later
	dt1 <- c(data,psum)
	ft <- eval(formula,dt1)

	dv <- as.numeric(ft$predictor) - as.numeric(ft$response)
	if(! is.null(vm)) dv <- rm %*% dv		
	sigma <- mad(dv)
	yold <- max(y)

									#///*************************************************
	while (iterate <= maxiter)	#///*********     Start Iteration    ****************
	{
		o <- order(y) # don't need a full sort, only smallest and two largest
		theta <- theta[o,] # so could be done more efficiently
		th <- th[o,]
		ht<-ht[o]
		y <- y[o]
		ilo <- 1
		ihi <- p+1
		inhi <- p
		rtol <- 2*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+1e-8)
		rtol2 <- sqrt(sum(th[1,]-th[2,])^2) / sqrt(sum(th[1,])^2)
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = y[ilo],th[ilo,],converge = min(rtol,rtol2)  ))
		if (rtol < tolerance | rtol2 < tolerance){
			Fault2 <- Fault()
			#result=list(parameters = as.list(th[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)
			break
		}
		if (iterate >= maxiter){
			Fault2 <- Fault(FN=1,FF = "optim.NM")
			#result=list(parameters = as.list(th[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)
			break
		}
		iterate <- iterate+2			# new point chosen by reflecting the worst current through the plane
											# of the others

		z <- smptry2(th,y,psum,data=data,ihi,-1,ht,vm=vm,rm=rm,sigma=sigma,control=control,robfunc=robfunc,formula=formula,...)
		if(is.Fault(z)) return(z)
		if (z[[1]] <= y[ilo]) {                                                         # new point is best--try going further
			z <- smptry2(z[[4]],z[[2]],z[[3]],data=data,ihi,2,ht,vm=vm,rm=rm,sigma=sigma,control=control,robfunc=robfunc,formula=formula,...)
			yold<-y[1]
			y <- z[[2]]; psum <- z[[3]]; th <- z[[4]];ht<-z[[5]]
			sigma <- z$sigma
		}
		else 
					if (z[[1]] >= y[inhi]) {
							ysave <- z[[2]][ihi]                                            #new point is still worst, try smaller step
							z <- smptry2(z[[4]],z[[2]],z[[3]],data=data,ihi,0.5,ht,vm=vm,rm=rm,sigma=sigma,control=control,robfunc=robfunc,formula=formula,...)
							yold<-y[1]
							y <- z[[2]]; psum <- z[[3]]; th <- z[[4]];ht<-z[[5]]            # a final psum is that after it not any smptry compute

							if (z[[1]] >= ysave) {                                          # still bad, shrink simplex
										for (i in (1:(p+1))[-ilo]) {
													psum <- (theta[i,]+theta[ilo,])/2
													theta[i,] <- psum
													th[i,] <- c(psum,sigma)
		
													if(is.null(vm)) 	ht[[i]] <- robloss(formula,data,start=th[i,] ,robfunc=robfunc,control=control,...)		## th: with sigma
													else ht[[i]] <- robloss.gn(formula,data,th[i,] ,robfunc,rm,control=control,...)
	
													if(is.Faultwarn(ht[[i]])){
														ht[[i]]$Fault$FN <- "optim.NL"
														return(ht[[i]])
													}
													y[i] <- as.numeric(ht[[i]]$htheta)
										}
										iterate <- iterate+p
										psum <- apply(th[,1:p],2,sum)
										dt1 <- c(data,psum)
										ft <- eval(formula,dt1)
										dv <- as.numeric(ft$predictor) - as.numeric(ft$response)
										if(! is.null(vm)) dv <- rm %*% dv		
										sigma <- mscale(dv)
							}
							else sigma <- z$sigma                                                # anytime the final psum must compute sigma
					} 
					else {
							y <- z[[2]]; psum <- z[[3]]; th <- z[[4]];ht<-z[[5]]
							sigma <- z$sigma
							iterate <- iterate-1
					}
	if(trace) lines(data[[formula$independent]],ht[[ilo]]$fmod$predictor)
	}	   	#### end itteration   #####
	##################################
#cat("finished iteration=",iterate ,"\n")
	#//////////*******************************************************************************************************************
	#//////////*******************************************************************************************************************
	#//////////********     preparing output   ***********************************************************************************
	#//////////********                        ***********************************************************************************
####	result=list(parameters = as.list(th[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)

	yresp <- as.numeric(ht[[ilo]]$fmod$response)
	ypred <- as.numeric(ht[[ilo]]$fmod$predictor)		
	ybar <- mean(yresp)	
	nlrho <- 1 - sum( ( yresp - ypred )  ^ 2 ) / 
					sum( (ypred - ybar ) ^ 2 )
#	curv1 <- curvature(gradient = attr(ht[[ilo]]$fmod$predictor,"gradient"),
#						hessian = attr(ht$fmod$predictor,"hessian"),
#						sigma = th[ilo,"sigma"])

if(is.null(vm))
		result <- nl.fitt.rob(
            parameters =   as.list(th[ilo,]),
		        scale = sigma,
						correlation =  nlrho,
						form =         formula,
						response =     ht[[ilo]]$fmod$response,
						predictor =    ht[[ilo]]$fmod$predictor,
#						curvature =    curv1,
						history =      iterhist,
						method = 		 fittmethod(methodID=			7,
														methodBR=			13,       ### rob or nopnrob variance
														detailBR=			"Nelder Mead",
														subroutine=		"nlmest.NM"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht[[ilo]]$htheta,
						rho =          ht[[ilo]]$rho,
						ri =           ht[[ilo]]$ri,
						curvrob =      NULL,
						robform =      robfunc
					)
	else

		result <- nl.fitt.rgn(
            parameters =  as.list(th[ilo,]),
		        scale = sigma,
            correlation =  nlrho,
						form =         formula,
						response =     ht[[ilo]]$fmod$response,
						predictor =    ht[[ilo]]$fmod$predictor, 
#						curvature =    curv1,
						history =      iterhist,
						method = 		 fittmethod(methodID=			8,
														methodBR=			13,       ### rob or nopnrob variance
														detailBR=			"Nelder Mead",
														subroutine=		"nlmest.NM"),
						data =         as.list(data),
						sourcefnc =     match.call(),
						Fault =        Fault2,
						htheta =       ht[[ilo]]$htheta,
						rho =          ht[[ilo]]$rho,
						ri =           ht[[ilo]]$ri,
						curvrob =      NULL,
						robform =       robfunc,
						vm =           vm,
						rm=            rm,
						gresponse =    transformNR(ht[[ilo]]$fmod$response,rm),
						gpredictor =   transformNR(ht[[ilo]]$fmod$predictor,rm)
					)
	return(result)
}
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of nlmest.NM'                              | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, Dep Stat, Stock                    | *
#* |                                                                        | *
#* |                Recoded from 'optim.NM' in March 20103                  | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************


	##
	##
	######  function smptry     ######
	
smptry2 <- function(th,y,psum,data,ihi,fac,ht,vm,rm,sigma,formula,control,...) {       #### sigma is the scale computed at psum
	ndim <- ncol(th)-1
	n <- length(data[[1]])

	fac1 <- (1-fac)/ndim
	fac2 <- fac1-fac
	ptry <- psum*fac1-th[ihi,1:ndim]*fac2
	th2 <- as.list(ptry)
	th2$sigma <- sigma
	if(is.null(vm)) htobj <- robloss(formula,data,start=th2 ,...)		## th: with sigma
	else htobj  <- robloss.gn(formula,data,th2,rmat=rm,control=control,...)
  if (is.Fault(htobj)) return(htobj)
	tmp <- as.numeric(htobj$ri)

	if(control$robscale) sigma<-mscale(tmp) # <- nl.mscale(tmp,robfunc,...)
	else sigma  <- sum(abs(tmp)) / n


	ytry <- as.numeric( htobj$htheta)
	if (ytry < y[ihi]) {
		y[ihi] <- ytry
		psum <- psum-th[ihi,1:ndim]+ptry
		th[ihi,1:ndim] <- ptry
		ht[[ihi]]<-htobj

		dt1 <- c(data,psum)
		ht0 <- eval(formula,dt1)
		dv <- as.numeric(ht0$predictor) - as.numeric(ht0$response)
		if(! is.null(vm)) dv <- rm %*% dv		
		sigma <- mscale(dv)
	}
	list(ytry,y,psum=psum,th=th,ht=ht,sigma=sigma)
}		
		
###############################################

