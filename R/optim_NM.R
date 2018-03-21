
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function optim.NLM, optimize function by Nelder Mead.                | *
#* |                                                                        | *
#* |    argumnts:                                                           | *
#* |      objfnc: any objective function for minimizing, it must contains:  | *
#* |              formula, data and start, extra will be defined in (...)   | *
#* |              the output of objfnc must contains:                       | *
#* |                $value(attr,gradient,hessian), $angmat,$angvec          | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      theta:   starting values, it must contains tau elements           | *
#* |      delta:   vector of size (p), include the name of parameters.      | *
#* |               for attaining (p+1) initial values.                      | *
#* |         Note: it depends on scale so its better to assigne by program  | *
#* |                                                                        | *
#* |                                                                        | *
#* |      ...:     can be entries objfnc                                    | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
optim.NM<-function(objfnc,data,start=getInitial(objfnc,data),delta=NULL,deltar=.1,
	control=nlr.control(tolerance=1e-8, maxiter=250 * length(start)),...)
{
	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	loc.start <- start
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	if(is.null(delta)){
		if(is.list(start)) delta<-unlist(start)*deltar
		else delta <- start*deltar
    delta[delta==0] <- deltar
	}
	else{
	  delta<-delta[names(start)]
	}
	p <- length(start)
	n <- length(data[[1]])
	th <- loc.start
	b<-unlist(start)
	thetaoptim <- rep(b,p)
	dim(thetaoptim) <- c(p,p)
	thetaoptim <- t(thetaoptim)
	diag(thetaoptim) <- diag(thetaoptim)+delta
	thetaoptim <- rbind(b-delta,thetaoptim)
  colnames(thetaoptim) <- names(start)
	y <- rep(0,p+1)
	iterhist <- NULL
	ht <- list(NULL)
  for (i in 1:(p+1)){
		ht[[i]] <- objfnc(data=data,start=thetaoptim[i,],...)
		if(is.Faultwarn(ht[[i]])){
			ht[[i]]@FF <- "optim.NM"
			return(ht[[i]])
		}
		y[i] <- as.numeric(ht[[i]]$value)
		iterhist <- rbind(iterhist,c(iteration = 0,objfnc = as.numeric(ht[[i]]$value),
				thetaoptim[i,],converge = 0))	
	}
	eol <- F
	iterate <- 0
	psum <- apply(thetaoptim,2,sum)
									#///*************************************************
	while (iterate <= maxiter)	#///*********     Start Iteration    ****************
	{
		o <- order(y) # don't need a full sort, only smallest and two largest
		thetaoptim <- thetaoptim[o,] # so could be done more efficiently
    if(p==1) {
			thetaoptim <- as.matrix(thetaoptim)
			colnames(thetaoptim) <- names(start)
		}
		ht<-ht[o]
		y <- y[o]
		ilo <- 1
		ihi <- p+1
		inhi <- p
		difference <- 2*abs(y[ilo+1]-y[ilo])
		rtol <- 2*abs(y[ilo+1]-y[ilo])/(abs(y[ilo])+1e-8)
  	iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = y[ilo],	thetaoptim[ilo,],converge = rtol))
    if (difference < tolerance*abs(y[ilo])){ 
			result=list(parameters = as.list(thetaoptim[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)
			return(result)
		}
		if (iterate >= maxiter){
			Fault2 <- Fault2<-Fault(FN=1,FF = "optim.NM")
			result=list(parameters = as.list(thetaoptim[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)
			return(result)
		}
		iterate <- iterate+2			# new point chosen by reflecting the worst current through the plane
										# of the others
		z <- smptry(thetaoptim=thetaoptim,y=y,psum=psum,objfnc=objfnc,data=data,ihi=ihi,fac=-1,ht=ht,...)
		if (z[[1]] <= y[ilo]) {		 # new point is best--try going further
			z <- smptry(thetaoptim=z[[4]],y=z[[2]],psum=z[[3]],objfnc=objfnc,data=data,ihi=ihi,fac=2,ht=ht,...)
			if(is.Fault(z)) return(z)
			y <- z[[2]]; psum <- z[[3]]; thetaoptim <- z[[4]];ht<-z[[5]]
		}
		else
			if (z[[1]] >= y[inhi]) {
				ysave <- z[[2]][ihi] #new point is still worst, try smaller step
				z <- smptry(thetaoptim=z[[4]],y=z[[2]],psum=z[[3]],objfnc=objfnc,data=data,ihi=ihi,fac=0.5,ht=ht,...)
				if(is.Fault(z)) return(z)
				y <- z[[2]]; psum <- z[[3]]; thetaoptim <- z[[4]];ht<-z[[5]]
				if (z[[1]] >= ysave) { # still bad, shrink simplex
					for (i in (1:(p+1))[-ilo]) {
						psum <- (thetaoptim[i,]+thetaoptim[ilo,])/2
						thetaoptim[i,] <- psum
						ht[[i]] <- objfnc(data=data,start=thetaoptim[i,],...)
						if(is.Faultwarn(ht[[i]])){
							ht[[i]]$Fault$FN <- "optim.NL"
							return(ht[[i]])
						}
						y[i] <- as.numeric(ht[[i]]$value)
					}
					iterate <- iterate+p
					psum <- apply(thetaoptim,2,sum)
				}
			} 
			else {
					y <- z[[2]]; psum <- z[[3]]; thetaoptim <- z[[4]];ht<-z[[5]]
					iterate <- iterate-1
			}
	}	   	#### end itteration   #####
	##################################
	## Creating output object
	##
	
	Fault2 <- Fault(FN=1,FF = "optim.NM")
  if(iterate>=maxiter) Fault2<-Fault(FN=1,FF = "optim.NM")
	iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = y[ilo],
										thetaoptim[ilo,],converge = -1)
						)
	result=list(parameters = as.list(thetaoptim[ilo,]), objfnc=ht[[ilo]], history=iterhist,Fault=Fault2)
	return(result)
}
	########################
		#*------------------------------------------*#
 		#*  end optim_NM                            *#
		#*  Hossein Riazoshams                      *#
		#*  15-Aug-2012                             *#
		#*  reference: statistical computing by R   *#
		#*------------------------------------------*#

	######  function smptry
smptry <- function(thetaoptim,y,psum,objfnc,data,ihi,fac,ht,...){
	ndim <- ncol(thetaoptim)
	fac1 <- (1-fac)/ndim
	fac2 <- fac1-fac
	ptry <- psum*fac1-thetaoptim[ihi,]*fac2
	thetaoptim2 <- thetaoptim
	y2 <- y
	ht2 <- ht
	psum2 <- psum
	htobj <- objfnc(data=data,start=ptry,...)
#	cat("psum,fac1,fac2,ndim,ptry \n",psum,fac1,fac2,ndim,ptry,"\n", 
#	"htobj$value   =",htobj$value,"\n")
	if(is.Fault(htobj)) return(htobj)
	ytry <- as.numeric( htobj$value)
	if (ytry < y2[ihi]) {
		y2[ihi] <- ytry
		psum2 <- psum-thetaoptim[ihi,]+ptry
		thetaoptim2[ihi,] <- ptry
		ht2[[ihi]]<-htobj
	}
	return(list(ytry,y2,psum2,thetaoptim2,ht2))
}
