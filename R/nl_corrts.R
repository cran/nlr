#******************************************************************************************************

nl.corrts <- function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.0010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),
    correlation=NULL,...)

{
	
	tols1 <- nlsqr(formula, data=data, start=start,control=control,...)
	if(is.Fault(tols1)) return(tols1)
	ri <- residuals(tols1) # should be corrected to work with formula in spherical norm
  n <- length(ri)
  switch(class(correlation)[1],
      "corAR1"={
        tm <- ar(ri,order.max=1,aic=F)
        cs2<-Initialize(correlation,data=as.matrix(ri))
        vmat <- corMatrix(cs2)
        #vinv <- solve(vmat)
        #umat <- chol(vmat)
        #ut <- t(umat)
        #rmat <- solve(ut)
        
        .temp1 <- 1.0/(1.0-tm$ar^2)
        vinv <- diag(c(1,rep(1+tm$ar^2,n)))
        vinv[col(vinv)==row(vinv)+1] <- -tm$ar
        vinv[row(vinv)==col(vinv)+1] <- -tm$ar
        rmat <- diag(c(sqrt(1-tm$ar^2),rep(1,n-1)))
        rmat[row(rmat)==col(rmat)+1] <- -tm$ar
        rmat <- rmat / sqrt(1-tm$ar^2)
        autpar<-tm$ar
      },
      
      "corARMA"={
        pcorr <- attr(correlation,"p")
        qcorr <- attr(correlation,"q")
        ncorr <- pcorr+qcorr
        tm <- arima(ri,order=c(pcorr,0,qcorr),include.mean = FALSE)
        correst <- corARMA(tm$coef,form=attr(correlation,"formula"),p=pcorr,q=qcorr,fixed=attr(correlation,"fixed"))
        cs2<-Initialize(correst,data=as.matrix(ri))
        vmat <- corMatrix(cs2)
        v2 <- eiginv(vmat,symmetric=T,stp=F)
        if(is.Fault(v2)) return(v2)
        for(i in 1:n)
          for(j in i:n)
            v2[i,j] <- v2[j,i]
        umat <- chol(v2)
        ut <- t(umat)
        rmat <- solve(ut)
        autpar<-tm$coef
      },
      "corCAR1"={
        
      },
      "corCompSymm"={
        
      },
      
      "corExp"={
        
      },
      "corGaus"={
        
      },
      "corLin"={
        
      },
      "corRatio"={
        
      },
      "corSpher"={
        
      },
      "corSymm"={
        
      }
      
    )
	tolerance <- control$tolerance*1e3
	minlanda<-control$minlanda / 1e4
	t2st <- nlsqr.gn(formula,data=data,
                   start=tols1$parameters[names(formula$par)],
                   vm=vmat, rm=rmat,
	                 control=nlr.control(tolerance=tolerance,minlanda=minlanda),
	                 ...)
	if(is.Fault(t2st)) return(t2st)
	t2st@autpar<-as.list(autpar)
	#t2st@autcorr <- tm
	return(t2st)
}
