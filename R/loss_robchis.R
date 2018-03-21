#************************************************************************
#**   hetroscedastic loss unction, minus of pseudo function            **
#**   this function compute the minus of pseoudo likelihood,based on   **
#**   model function.                                                  **
#**   parameters and data vlues must be stored as bellow.              **
#**   to find l(tau)=sum[log(h(mu;tau,sg))] + sum[rho(ri/h(mu;tau,sg))]**
#**     mu = y~f(x,^theta)                                             **
#**     ri = y-f(x,^theta)                                             **
#**   parameters and and ariables.                                     **
#**    Data.                                                           **
#**     (xr,yr,t) stored in data list, as general can be stored.       **
#**     Note:   data must be ordered.                                  **
#**     t=mu, if is null must me computed, otherwise the user entry co **
#**    theta.   theta list.                                            **
#**    tau.                                                            **
#**     (sg,landa) or any other form organized by user.                **
#**    the varmodel must return back                                   **
#**            stdev = h = sigma*g ------standard dev.--               **
#**            variance will be change to sigma^2 g^2                  **
#**                2 h grad(h)                                         **
#**                2 h gr(h)T grd(h) + 2 h hes(h)                      **
#**    Note: The entry of variance is (f) not (x), H(f)                **
#************************************************************************
loss.robchis <- function(formula,data,start,theta,varmodel,robfunc,...){
	Fault2 <- Fault()
	data <- as.list(data)
	theta <- as.list(theta)
	datalist <- c(data,theta)
	fmod <- eval(formula,datalist)                  #  f(x) n*1                           #
	pred <- as.numeric(fmod$predictor)              #  f(x) values                        #
	resp <- as.numeric(fmod$response)               #  y, left side                       #
	ri <- resp - pred
	nrp <- nonrepl(list(x=data$xr,y=data$yr))
	data$xr <- nrp$x
	data$yr <- nrp$y
	z <- rzvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]
	zi=z
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
	if(is.Fault(varcomp)) return(Fault(FN=17))
	vc <- as.numeric(varcomp$predictor)
	vch <- attr(varcomp$predictor,"hessian") 
	vcg <- attr(varcomp$predictor,"gradient")
	z <- z[nrp$xm]
	vc <- vc[nrp$xm]
	vch <- vch[nrp$xm,,]
	vcg <- vcg[nrp$xm,]
	wi <- pmax(1,nrp$ni-1)
	refvar <- sum(wi * z)/sum(wi)
	vc2<-vc
#	vc2 <-pmax(vc,0.01*refvar)
	vc[vc<0] <- 0.05 #0.1*refvar
	s1 <- wi * log(vc2)                              #  (ni-1) log(g)                     #
	s2 <- sqrt(z / vc2)                              #   sqrt(zi) / gi (gi standard devia)#
	robvalue <- robfunc$fnc(s2,...)                 #   rho(zi/gi)
	rbv <- as.numeric(robvalue)
	rbg <- attr(robvalue,"gradient")
	rbh <- attr(robvalue,"hessian")
	value <- sum(s1+wi*rbv)                         #  sum((ni-1)log(gi)+wirho(\/zi/\gi) #
	temp2 <- 0.5 * wi * sqrt(z/vc^3)                #   zi / gi^2                        #
	temp3 <- wi / vc                        #   (ni-1) / gi                      #
	temp5 <- temp3 - temp2 * rbg                    #  2(ni-1)/gi - zi/gi^2.rh           #
	angvec <- temp5
	temp4 <- t(vcg) %*% temp5                       #  gradg * [2(ni-1)/gi - zi/gi^2.rh ]#
	gradient <- temp4                               #   gradinet (like)                  #
	temp4 <- prodVA(vch,temp5)                      #  hess(g) * 2(ni-1)/gi - zi/gi^2    #
	.expr1 <- (vc^5.0)
	.expr2 <- wi * z
	.expr3 <- ( .expr2/ sqrt(.expr1) )
	temp6 <- .expr3 * rbg * (3.0/4.0)               #  (  3 * wi*zi/ 4*sqrt(gi^5))*rho.  #
	temp7 <- wi / vc^2                              #    (ni-1) / gi^2                   #
	exm1 <- (wi*z) / (4*vc^3)                       # (wi zi/ 4.gi^3)                    # 
	exm2 <- exm1 * rbh                              # (wi zi/ 4.gi^3) * rho(zi/gi)       # 
	a <- temp6-temp7
	temp8 <- a +exm2										# 2 * zi/gi^3  - (ni-1) / gi^2       #
	temp9 <- diag(temp8)                            # diag(2.zi/gi^3  - 2.(ni-1)/gi^2+(z/g2)2rh)   #

	t1 <- t(vcg) %*% temp9 	                        # grad * diag(2 * zi/gi^3  - (ni-1) / gi^2)    #

	t2 <- t1 %*% vcg                                # grad*diag(2 * zi/gi^3  - (ni-1) / gi^2+..rh)*grd  #
	hessian <- temp4 + t2                           #  hessian (likelihood)                        #
 	gradient=t(gradient)
	vcmdata <- list()
	vcmdata[[varmodel$independent]]<-data[[varmodel$independent]]
	vcmdata[[varmodel$dependent]]<-data[[varmodel$dependent]]	
	attr(value,"gradient") <- gradient
	attr(value,"hessian") <- hessian	
	result <- list (value=value,angvec=angvec,angmat=vcg,
		refvar=refvar,fmod=fmod,varcomp=varcomp,vcmdata=vcmdata,sourcefnc= match.call(),
		rho=robvalue,zi=zi)
	return(result)
	#****************************************************
	#*****   output:                                ***** 
	#*****  likelihood of chisquare for variance.   ***** 
	#*****   and its gradient and hessian.          *****
	#*****  standard output for optim.LM or others  *****
	#*****                                          *****
	#*****     writhen by Hossein Riazoshams        *****
	#*****              08/01/2010                  *****
	#*****                                          *****
	#****************************************************

}

