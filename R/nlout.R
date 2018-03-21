
#+------------------------------------------------------------------------------+
#|    nlout: find the measurments for outliers.                                 |
#|       nlfited: nl.fitt or nl.fitt.rob or nl.fitt.mm or any inherited object. |
#|    detail. the residuals and x can be robust but it make confusion.          |
#+------------------------------------------------------------------------------+
nlout<-function(nlfited){

  if(class(nlfited)=="Fault") return(nlfited)
	if(nlfited$Fault$FN!=0) return(nlfited$Fault)
	v <- attr(nlfited@predictor,"gradient")
	hs <- attr(nlfited@predictor,"hessian")
	res <- residuals(nlfited)
	sigma <- nlfited@parameters[["sigma"]]
	grtgr <- t(v) %*% v									##   F.' * F.     P*P
	ht<-eiginv(grtgr,symmetric =T)						##   (F.' F.)^-1  P*P
	vbar <- apply(v,2,mean)								##   Mean over columns of v, 
	if(class(ht)=="Fault") return(ht)
	covmat <- var(v)										##   (C) Variance covariance matrix.
	cinv <- eiginv(covmat,symmetric =T)
	p <- nlfited$form$p
	n <- length(res)
	wnew <- potmah <- vmat <- mahd <- rep(0,n)

	dataset <- nlfited$data[c(nlfited$form$independent,nlfited$form$dependent)]	
	dataset <- data.frame(dataset)
	g2<-cov.mve(dataset)
	mahdata <- mahalanobis(dataset,center=g2$center,cov=g2$cov)
	mahdata <- sqrt(mahdata)
	ctfmahd1 <- median(mahdata) + 3 * mad(mahdata)	##  Mahalanobis and ctfp for all data together

	R <- (mahdata <= ctfmahd1)
	xr <- as.matrix(dataset[R,])						##  X, whole data, Remaining
	xtx <- t(xr) %*% xr									##  XT X
	wninv <- eiginv(xtx)									## (XT X) ^-1
	datasetx <- nlfited$data[c(nlfited$form$independent)]
	datasetx <- data.frame(datasetx)
	g3<-cov.mve(datasetx)
	mahdatax <- mahalanobis(datasetx,center=g3$center,cov=g3$cov)
	mahdatax <- sqrt(mahdatax)
	ctfmahd2 <- median(mahdatax) + 3 * mad(mahdatax) ##  Mahalanobis for only x data

															##  compute vmat: square of mahalanobis distance
	for(i in 1:nrow(v)){									##   Cllasic, center and variance are classic.
		b1<-v[i,]											##  vi', will store in row, 1*p
		b2<-b1 %*% ht										##  vi' (g'g)^-1            1*p
		vmat[i] <- b2 %*% b1								##  wii = vi' (g'g)^-1 vi   1*1
		temp <- v[i,] - vbar								##  vi - mean(v)            1*p
		temp2 <- temp %*% cinv							##  (vi - vbar) * C^-1      p*1, first is row
		temp2 <- temp2 %*% (temp)						##  (vi-vbar) C^-1 (vi-vbar) 1*1 second column, 
															##                                first row
		mahd[i] <- sqrt(temp2)							##  result after sum is a vectors
		b1 <- as.matrix(dataset[i,])
		b2<- b1 %*% wninv									##  vi' (xr'xr)^-1            1*p
		wnew[i] <- b2 %*% t(b1)							##  vi' (xr'xr)^-1 vi         1*1
	}
	ctfmahdr <- median(mahd)+3*mad(mahd)
	studres <- res / sigma / sqrt(1-vmat)				##  studentized residual, ri / (sigm sqrt(1-wii))
	elliptnorm <- studres^2 * vmat / (1-vmat) / p	##  Influence curve, Elliptic norm (Cook d)
	hadi <- vmat / (1-vmat)								##  Hadi to assess high leverage points wii/(1-wii)
	mdhad <- median(hadi)
	ctfhadi1 <- mdhad + 2 * mad(hadi)
	ctfhadi2 <- mdhad + 3 * mad(hadi)
	potmah[R] <- wnew[R] / (1 - wnew[R])				##  Potential maha lanobis Remain = wii(-d)/(1-wii(-d))
	potmah[!R] <- wnew[!R]								##  Potential for deleted data

	ctfpotmah <- median(potmah) + 3 * mad(potmah)	##  

															#####################################
	JL <- JacobianLeverage(nlfited)					###   Jacobian Leverage
	if(class(JL)!="Fault"){
		jvmat <- diag(JL)
		jl.studres <- res / sigma / sqrt(1-jvmat)					##  studentized residual, ri / (sigm sqrt(1-wii))
		jl.elliptnorm <- jl.studres^2 * jvmat / (1-jvmat) / p	##  Influence curve, Elliptic norm (Cook d)
		jl.hadi <- jvmat / (1-jvmat)								##  Hadi to assess high leverage points wii/(1-wii)
		mdhad <- median(jl.hadi)
		jl.ctfhadi1 <- mdhad + 2 * mad(jl.hadi)
		jl.ctfhadi2 <- mdhad + 3 * mad(jl.hadi)
	}
	###############################################
   ##      Single deletion                      ##
	###############################################
	d.stud <- d.sigma <- d.yhat <- d.ffits <- rep(0,n)
	d.fbetas <- matrix(rep(0,n*p),nrow=n)
	yi <- as.numeric(nlfited$predictor)
	yihat <- as.numeric(nlfited$response)
			############################################### new idea for deleted object   #############

	for(i in 1:n){
		dst2 <- dataset[-i,]									##   xy(-i)
		newfiti <- recalc(nlfited,dst2)					##   fit all (-i)
		newfiti@Fault@FF <- "nlout"
		if(newfiti$Fault$FL) return(newfiti)
		d.sigma[i]  <- newfiti$parameters[["sigma"]] 	##  sigma(-i)
		yhati <- predict(newfiti,newdata=dataset[i,])				##  yhat(-i) grd&hess
		d.yhat[i] <- as.numeric(yhati)						##  yhat(-i)
		dfbt <- unlist(nlfited$parameters[names(nlfited$parameters)!="sigma"]) - 
		        unlist(newfiti$parameters[names(newfiti$parameters)!="sigma"])
		d.fbetas[i,] <- dfbt / newfiti$parameters[["sigma"]] / sqrt(diag(ht))   ##   DFBETASi = bj-bj(i) / sg(i) sart(x'x-1 jj)
	}
	delstud <- (yi - d.yhat) / d.sigma / sqrt(1-vmat)					##  deleted resid
	yhdiff <- yihat - d.yhat												##  yi^ - y(-i)^
	d.ffits <- abs(delstud) * sqrt( vmat / (1-vmat) )					##  dffits
	ctfdffits <- 2 * sqrt(p/n)												##  cut of poin dffits
	atkinson <- sqrt( (n-p) / p ) * abs(d.ffits)  						##  Atkinsons distance
	ctfatk <- median(atkinson) + mad(atkinson)
	if(class(JL)!="Fault"){
		jl.delstud <- (yi - d.yhat) / d.sigma / sqrt(1-jvmat)				##  deleted resid
		yhdiff <- yihat - d.yhat												##  yi^ - y(-i)^
		jl.d.ffits <- abs(delstud) * sqrt( jvmat / (1-jvmat) )				##  dffits
		ctfdffits <- 2 * sqrt(p/n)												##  cut of poin dffits
		jl.atkinson <- sqrt( (n-p) / p ) * abs(jl.d.ffits) 					##  Atkinsons distance
		jl.ctfatk <- median(jl.atkinson) + mad(jl.atkinson)
	
}
	################################################################### output   #########################
	result1 <- structure(.Data=list(
								vmat,
								d.yhat,
								nl.robmeas(measure=studres,cutofpoint=c(2.5,3,-2.5,-3),name="Studentised Residuals"),
								nl.robmeas(measure=sqrt(elliptnorm),cutofpoint=1,name="Elliptic Norm (Cook Dist)"),
								nl.robmeas(measure=mahd,cutofpoint=c(ctfmahdr,qchisq(.95,p)),name="Regression Mahalanobis Distance"),
								nl.robmeas(measure=mahdata,cutofpoint=ctfmahd1,name="Mahalanobis MVE, data"),
								nl.robmeas(measure=mahdatax,cutofpoint=ctfmahd2,name="Mahalanobis MVE, xs"),
								nl.robmeas(measure=hadi,cutofpoint=c(ctfhadi1,ctfhadi2),name="Hadi potential"),
								nl.robmeas(measure=potmah ,cutofpoint=ctfpotmah,name="Potential mahalanobis"),
								nl.robmeas(measure=delstud ,cutofpoint=c(2.5,3,-2.5,-3),name="Deletion Studentized"),
								nl.robmeas(measure=d.ffits ,cutofpoint=ctfdffits ,name="DFFITS"),
								nl.robmeas(measure=atkinson ,cutofpoint=c(2,ctfatk),name="Atkinson Distance"),
								g2,g3,
								d.fbetas
								
							  ),
							.Names=c(
										"vmat",
										"yihat",
										"studres",
										"cook",
										"mahd.v",
										"mahd.dt",
										"mahd.xs",
										"hadi",
										"potmah",
										"delstud",
										"dffits",
										"atk",
										"mvedta","mvex",
										"dfbetas")
							)

	if(class(JL)!="Fault"){ 
		result2 <- structure(.Data=list(
								jvmat,
								nl.robmeas(measure=jl.studres,cutofpoint=c(2.5,3,-2.5,-3),name="Jacobian Leverage Studentised Residuals"),
								nl.robmeas(measure=sqrt(jl.elliptnorm),cutofpoint=1,name="Jacobian Leverage Elliptic Norm (Cook Dist)"),
								nl.robmeas(measure=jl.hadi,cutofpoint=c(jl.ctfhadi1,jl.ctfhadi2),name="Jacobian Leverage Hadi potential"),
								nl.robmeas(measure=jl.delstud ,cutofpoint=c(2.5,3,-2.5,-3),name="Jacobian Leverage Deletion Studentized"),
								nl.robmeas(measure=jl.d.ffits ,cutofpoint=ctfdffits ,name="Jacobian Leverage DFFITSi"),
								nl.robmeas(measure=jl.atkinson ,cutofpoint=c(2,jl.ctfatk),name="Jacobian Leverage Atkinson Distance")
								),
							.Names=c(
										"jl.vmat",
										"jl.studres",
										"jl.cook",
										"jl.hadi",
										"jl.delstud",
										"jl.dffits",
										"jl.atk")
							)
		result1[names(result2)] <- result2
	}
	
	return(result1)
}

#--------------------------------------------------------------------------------
#|                                  End of nlout. 
#|                                  29/12/2008
#--------------------------------------------------------------------------------

