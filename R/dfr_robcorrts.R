#******************************************************************************************************
dfr.robcorrts <- function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.0010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),correlation=1,
		robfunc,...)

{
	q <- correlation

	#	tols2 <- nl.mm3th(formula, data=data, start=start, nl.robfuncs[[2]])	
	tols2 <- nlmest.NM(formula, data=data, start=start,robfunc=robfunc,control=control,...)	
	if(is.Fault(tols2)) return(list(fited=tols2))
	ri <- residuals(tols2)
	st <- tols2$parameters

	tm <- ar(ri,order.max=q,aic=F,...)
#	tm <- ar.gm(ri,order=q,b=F)
	n <- length(data$xr)
	if(q==1){
		###### ----------  AR(1)
		.temp1 <- 1.0/(1.0-tm$ar^2)
		vinv <- diag(c(1, rep(1+tm$ar^2,n-1) ))
		vinv[col(vinv)==row(vinv)+1] <- -tm$ar

		vinv[row(vinv)==col(vinv)+1] <- -tm$ar
		vinv <- vinv / (1.0 - tm$ar^2)
		rmat <- diag(c(sqrt(1-tm$ar^2),rep(1,n-1)))
		rmat[row(rmat)==col(rmat)+1] <- -tm$ar
		rmat <- rmat / sqrt(1-tm$ar^2)
	}
	else if(q==2) {
		rho <- 1

		rho[2] <- tm$ar[1]/(1-tm$ar[2])
		for( i in 3:n)	rho[i] <- tm$ar[1] * rho[i-1] + tm$ar[2] * rho[i-2]
		vmat <- matrix(rep(0,n*n),nrow=n)
		for( i in 1:(n-1)){
			rho2 <-0
			rho2 <- rho[i:1]
			rho2[(i+1):n] <- rho[2:(n-i+1)]
			vmat[i,] <- rho2
		}
		vmat[n,] <- rho[n:1]
		v2 <- eiginv(vmat,symmetric=T,stp=F)
		if(is.Fault(v2)) return(list(fited=v2))
		for(i in 1:n)
			for(j in i:n)
				v2[i,j] <- v2[j,i]
		rmat <- chol(v2)

	}
	else if(q==3) {
			###### ----------  AR(3)
		rho <- 1

		rho[2] <- (tm$ar[1] + tm$ar[3] * tm$ar[2]) / ( 1-tm$ar[2] - tm$ar[3] *(tm$ar[1] + tm$ar[3] ) )
		rho[3] <- tm$ar[2] + (tm$ar[1] + tm$ar[3]) * rho[2]
		for( i in 4:n)	rho[i] <- tm$ar[1] * rho[i-1] + tm$ar[2] * rho[i-2] + tm$ar[3] * rho[i-3]
		vmat <- matrix(rep(0,n*n),nrow=n)
		for( i in 1:(n-1)){
			rho2 <-0
			rho2 <- rho[i:1]
			if(i!=n)rho2[(i+1):n] <- rho[2:(n-i+1)]
			vmat[i,] <- rho2
		}
		vmat[n,] <- rho[n:1]
		v2 <- eiginv(vmat,symmetric=T,stp=F)
		if(is.Fault(v2)) return(list(fited=v2))
		for(i in 1:n)
			for(j in i:n)
				v2[i,j] <- v2[j,i]
		rmat <- chol(v2)
	}
	t2st <- nlmest.NM(formula, data=data, start=tols2$parameters,robfunc= robfunc,vm=vmat,rm=rmat,control=control,...)	
	result <- list(fited=t2st,tm=tm)
}