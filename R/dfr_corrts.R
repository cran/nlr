#******************************************************************************************************

dfr.corrts <- function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.0010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),correlation=1,...)

{
	q <- correlation
	tols1 <- nlsnm(formula, data=data, start=start,control=control,...)
	if(is.Fault(tols1)) return(list(fited=tols1))
	ri <- residuals(tols1)
	tm <- ar(ri,order.max=q,aic=F,...)
	tm$ar<-as.numeric(tm$ar)
	n <- length(data$xr)
	if(q==1){
		###### ----------  AR(1)
		.temp1 <- 1.0/(1.0-tm$ar^2)
		vinv <- diag(c(1,rep(1+tm$ar^2,n)))
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
		v2 <- eiginv(vmat,symmetric=T)
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

	t2st <- nlsnm(formula,data=data,start=tols1$parameters,vm=vmat ,rm=rmat,...)
	return(list(fited=t2st,tm=tm))
}
