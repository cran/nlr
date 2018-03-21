#************************************************************************
#**   hetroscedastic loss unction, minus of pseudo function            **
#**   this function compute the minus of pseoudo likelihood,based on   **
#**   model function.                                                  **
#**   parameters and data vlues must be stored as bellow.              **
#**   to find l(tau)=sum[log(h(mu;tau,sg))] + sum[rho(ri/h(mu;tau,sg))]**
#**     mu = y~f(x,^thet)                                              **
#**     ri = y-f(x,^theta)                                             **
#**   parameters and and ariables.                                     **
#**    Data.                                                           **
#**     (xr,yr,t) stored in data list, as general can be stored.       **
#**     Note:   data must be ordered.                                  **
#**     t=mu, if is null must me computed, otherwise the user entry co **
#**    theta.   theta list.                                            **
#**    tau.                                                            **
#**     (sg,landa) or any other form organized by user.                **
#**    the varmodel must return back sigma^2*g --------variance.------ **
#************************************************************************
loss.chis <- function(formula,data,start,theta,varmodel,...){
	Fault2 <- Fault()
	data <- as.list(data)
	theta <- as.list(theta)
	datalist <- c(data,theta)
	fmod <- eval(formula,datalist)                  #  f(x) n*1                           #
	pred <- as.numeric(fmod$predictor)              #  f(x) values                        #
	resp <- as.numeric(fmod$response)               #  y, left side                       #
	ri <- resp - pred
	nrp <- nonrepl(list(x=data[[formula$independent]],y=data[[formula$dependent]]))
	data[[formula$independent]] <- nrp$x
	data[[formula$dependent]] <- nrp$y
	z <- zvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]        # variance z=zi , si^2                #
	zi=z
	######z <- sqrt(z)                              # standard deviation ????             #
	n <- length(resp)
	p <- varmodel$p
	if(is.null(data[[varmodel$independent]])) 
			data[[varmodel$independent]] <- pred     #  t=mu                               #
	data[[varmodel$dependent]] <- zi
	theta <- as.list(theta)
	tau <- as.list(start)
	datalist <- c(data,theta,tau)                   #  (x,y,t,theta,sg,landa)             #
	datalist=as.list(datalist)
	varcomp <- eval(varmodel,datalist)              #   H(sg,landa,t)                     #
	if(is.Fault(varcomp)) return(varcomp)
	vc <- as.numeric(varcomp$predictor)
	vch <- attr(varcomp$predictor,"hessian")        #   n*q*q                             #
	vcg <- attr(varcomp$predictor,"gradient")       #   n*q                               #
	z <- z[nrp$xm]
	vc <- vc[nrp$xm]	
	vch <- vch[nrp$xm,,]
	vcg <- vcg[nrp$xm,]
	wi <- pmax(1,nrp$ni-1)
	refvar <- sum(wi * z)/sum(wi)
	vc2=vc
#	vc2 <- pmax(vc,0.01*refvar)
	wi3 <- pmax(1,nrp$ni-3)
	s1 <- wi * log(vc2)                              #  (ni-1) log(g)                     #
#	s1 <- wi3 * log(vc2)                              #  (ni-1) log(g)                     #
	s2 <- z / vc2                                    #       (zi) / gi (gi standard devia)#
	value <- sum(s1+wi*s2)                          #   sum( (ni-1) log(gi)+zi/gi)       #
###	value <- sum(s1+s2)                             #   sum( (ni-1) log(gi)+zi/gi)       #
	refvar <- sum(wi * z)/sum(wi)
	temp2 <- (z / vc^2) * wi                        #   zi / gi^2                        #
	temp3 <- wi / vc                                #   (ni-1) / gi                      #
	temp5 <- temp3 - temp2                          #   (ni-1)/gi - zi/gi^2.rh           #
	angvec <- temp5
	temp4 <- t(vcg) %*% temp5                       #  gradg * [2(ni-1)/gi - zi/gi^2.rh ]#
	gradient <- temp4                               #   gradinet (like)                  #
	temp4 <- prodVA(vch,temp5)                      #  hess(g) * 2(ni-1)/gi - zi/gi^2    #
	temp6 <- ((2.0 * z) / vc^3) * wi                #      2.zi/gi^3                     #
###	temp7 <- (2.0 * wi) / vc^2                      #  2 (ni-1) / gi^2                   #
	temp7 <- wi / vc^2                              #  2 (ni-1) / gi^2                   #
	temp8 <- temp6 - temp7                          #     zi/gi^3  - (ni-1) / gi^2       #
	temp9 <- diag(temp8)                            # diag(zi/gi^3  - (ni-1)/gi^2        #
	t1 <- t(vcg) %*% temp9 	                        # grad * diag(zi/gi^3  - (ni-1) / gi^2)    #
	t2 <- t1 %*% vcg                                # grad*diag(zi/gi^3  - (ni-1) / gi^2+..rh) #
	hessian <- temp4 + t2                           #  hessian (likelihood)                        #
 	gradient=t(gradient)
	vcmdata <- list()
	vcmdata[[varmodel$independent]]<-data[[varmodel$independent]]
	vcmdata[[varmodel$dependent]]<-data[[varmodel$dependent]]
	attr(value,"gradient") <- gradient
	attr(value,"hessian") <- hessian	
	result <- list (value=value,angvec=angvec,angmat=vcg,
		refvar=refvar,fmod=fmod,varcomp=varcomp,vcmdata=vcmdata,sourcefnc= match.call(),zi=zi)
	return(result)


	#****************************************************
	#*****   output:                                ***** 
	#*****  likelihood of chisquare for variance.   ***** 
   #*****   and its gradient and hessian.          *****
   #*****  standard output for optim.LM or others  *****
   #*****                                          *****
   #*****     writhen by Hossein Riazoshams        *****
   #*****              22/01/2010                  *****
   #*****                                          *****
   #****************************************************

}





