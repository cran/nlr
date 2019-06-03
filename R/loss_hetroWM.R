#************************************************************************
#**   hetroscedastic loss unction, minus of pseudo function            **
#**   this function compute the minus of pseoudo likelihood,based on   **
#**   model function.                                                  **
#**   parameters and data vlues must be stored as bellow.              **
#**   to find l(tau,theta)=sum[log(h(mu;tau,sg))] +                    **
#**                                           sum[rho(ri/h(mu;tau,sg))]**
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
#**    Note: The entry of variance is (x) , H(x)                       **
#**          variance function is general function of xi.              **
#**          if want to be function of (f) it should be inside "data"  **
#************************************************************************
loss.hetroWM <- function(formula,data,start,varmodel,robfunc,...){
	Fault2 <- Fault()
	data <- as.list(data)
	start<-as.list(start)

  theta <- start[names(formula$par)]
	theta <- as.list(theta)
	tau <- start[names(varmodel$par)]
	tau<-as.list(tau)
	datalist	 <- c(data,start)                     #  (x,y,t,theta,sg,landa)             #
	datalist=as.list(datalist)
	fmod <- eval(formula,datalist)                  #  f(x) n*1                           #
	pred <- as.numeric(fmod$predictor)              #  f(x) values                        #
	predg <- attr(fmod$predictor,"gradient")        #  f'   gradient                      #
	predh <- attr(fmod$predictor,"hessian")         #  f"   hessian                       #
	resp <- as.numeric(fmod$response)               #  y, left side                       #
	ri <- resp - pred
	n <- length(resp)
	q <- varmodel$p
	p <- formula$p
	sigma <- mad(ri)# sum(ri^2) / (n-p)
	datalist[[varmodel$dependent]] <- ri

	if(is.null(data[[varmodel$independent]])){
		datalist[[varmodel$independent]] <- pred        #  t=mu                           #
	} 
	varcomp <- eval(varmodel,datalist)              #   H(sg,landa,t)                    #
	if(is.Fault(varcomp)) {
		return(Fault(FN=17,FF="loss.hetroWM"))
	}
	vc <- as.numeric(varcomp$predictor)
	vcg <- attr(varcomp$predictor,"gradient")
	vch <- attr(varcomp$predictor,"hessian")
	nrp <- nonrepl(list(x=data$xr,y=data$yr))
	z <- rzvalues(ri,nrp$ni,nrp$xo) #[nrp$xm]
	zi=z
	z <- z[nrp$xm]
	wi <- pmax(1,nrp$ni-1)
	refvar <- sum(wi * z)/sum(wi)
	vc[vc<0] <- 0.05# 0.01*refvar
	
	#********* calculate standard deviation and all derivatives          ****************#
	#********* This is transforming from variance to standard deviation  ****************#

	sd <- sqrt(vc)                                  #  sigma i (standard deviation)      #
	sdg <- -vcg / (2*sd)                            #  grad (sigmai)   gradient stdev    #
	vcgt <- t(vcg)
	gtg <- vcgt %m3d% vcg                           # vcgT * vc (n*p*p)                  #
	sd3 <- sd * vc
	.temp2 <- 4.0 * sd3                             # 4* sigma ^3                        #
	.expr2 <- (1.0 / .temp2) * gtg                  # 1/4* sigma^3 * gr T * grad         #
	.expr1 <- (0.5 / sd) * vch                      # 1/2*sigma  * hess                  #
	sdh <- -.expr1 + .expr2

   #************************   the value ***********************************************#

	s1 <- log(sd)                                   #  (ni-1) log(g), instead used log(g)#
	s2 <- (ri) / sd                                 #   ri / gi (gi var)# 
	robvalue <- robfunc$fnc(s2,...)                 #   rho(ri/gi)
	rbv <- as.numeric(robvalue)
	rbg <- attr(robvalue,"gradient")
	rbh <- attr(robvalue,"hessian")
	value <- sum(rbv+s1)                             #  sum((ni-1)log(gi)+wirho(\/zi/\gi)#

   # ****************************         gradient     ************************************#
	.midle1 <- rbg / vc                                #  psi / hi                          #
	.minmidle <- -.midle1                              #  -psi / sigmai                     #
	grad.theta <- t(predg) %*% .minmidle               #  g(f)T * (-psi / hi)               #
	rbgri <- rbg * ri
	.temp3 <- rbgri / vc                            #  -psi / sigma^2                    #
	.temp4 <- (1.0 / sd) - .temp3                   #
	grad.tau <- t(sdg) %*% .temp4
	gradient <- rbind(grad.theta,grad.tau)        # gradinet, both theta & tau         #
	gradient <- t(gradient)
	               #******    convergance criterion *************************************#
	
	.z1 <- matrix(rep(0,n*p),ncol=p)
	.z2 <- matrix(rep(0,n*q),ncol=q)
	angmat <- rbind(cbind(predg,.z2),cbind(.z1,sdg))
	angvec <- c(.temp2,.temp4)

                  #********   hessian              *************************************#
	
	.midle2 <- rbh / vc^2
	.temp2 <- diag(.midle2)
	.temp3 <- t(predg) %*% .temp2
	.expr1 <- .temp3 %*% predg
	.expr2 <- prodVA(predh,.midle1)                  #  hess(f) * h'/sigma                #
	hrow11 <- .expr1 - .expr2
	
	.temp1 <- .midle2 + ( 2 * (rbg/sd3) )
	.temp2 <- (t(predg)) %*% (diag(.temp1))
	hrow12 <- .temp2 %*% sdg                         # g(f)T * [h"*ri/sigma^3 + h'/sigma^2#
	hrow21 <- t(hrow12)
	
	.t1 <- rbh * (ri^2)
	.temp1 <- .t1  / (vc * vc)
	.t2 <- rbg * ri
	.temp2 <- (2 * .t2) / sd3
	.temp3 <- 1.0 / vc
	.midle3 <- .temp1 + .temp2 - .temp3
	.temp4 <- diag(.midle3)
	.temp4 <- (t(sdg)) %*% .temp4
	.expr1 <- .temp4 %*% sdg

	.temp1 <- 1.0/sd
	.t1 <- rbg * ri
	.t2<- .t1 / vc
	.midle4 <- .temp1 - .t2
	.expr2 <- prodVA(sdh,.midle4)
	hrow22 <- .expr1 + .expr2
	hessian <- rbind(cbind(hrow11,hrow12),cbind(hrow21,hrow22)) #  hessian (likelihood)  #
	
	attr(value,"gradient") <- gradient
	attr(value,"hessian") <- hessian	
	nlrho <- 1 - sum( ( resp - pred )  ^ 2 ) / 
					sum( (pred - mean(resp) ) ^ 2 )

	result <- list (value=value,angvec=angvec,angmat=angmat,
		refvar=refvar,sourcefnc= match.call(),
			rho=robvalue,fmod=fmod,varcomp=varcomp,correlation =nlrho,ri=ri)
	return(result)
	#****************************************************
	#*****   output:                                ***** 
	#*****  likelihood of rob M (normal)            ***** 
	#*****   and its gradient and hessian.          *****
	#*****  standard output for optim.LM or others  *****
	#*****                                          *****
	#*****     writhen by Hossein Riazoshams        *****
	#*****              02/07/2012                  *****
	#*****                                          *****
	#****************************************************

}

