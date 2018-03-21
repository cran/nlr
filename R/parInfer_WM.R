parInfer.WM <- function(object,confidence=.95){

	sdcomp <- transfquad(object$hetro$predictor)
	vvar <- as.numeric(object$hetro$predictor)
	vsdv <- as.numeric(sdcomp)
	ka <- 1.0 / vsdv
	gvar <- attr(sdcomp,"gradient")
	hvar <- attr(sdcomp,"hessian")
	epsilon <- residuals(object) / vsdv
	
	rhotcvt <- (object$rho)
	
	psi <- attr(rhotcvt,"gradient")
	psid <- attr(rhotcvt,"hessian")
	sgpsi2 <- mean(psi^2)
	vrpsi2 <- var(epsilon * psi)
			
	gama1 <- mean(epsilon * psi)
	gama2 <- mean(psid)
	.temp <- epsilon * epsilon
	gama3 <- mean(.temp * psid)
	
	gradf <- attr(object$predictor,"gradient")
	n <- nrow(gradf)
	p <- object$form$p
	q <- object$hetro$form$p
	
	fft <-  t(gradf) %m3d% gradf                            ### f * fT  n*p*p
	ka2 <- ka^2

	kfft <- prodVA(fft,ka2)
	G1n <- gama2 * kfft

	G31n <- sgpsi2 * kfft
	ss <- t(gvar) %m3d% (gvar)

	G2n <- matrix(rep(0,q*q),nrow=q)
	for (i in 1:n){
		ad1 <- ((2*gama1 + gama3 -1) / vvar[i]) * ss[i,,]
		ad2 <- ((1-gama1)/vsdv[i]) * hvar[i,,]
		ad <- ad1+ad2
		G2n <- G2n + ad
	}
	kvvt <- prodVA(ss,ka2)
	G32n <- vrpsi2 * kvvt
	
	zero1 <- matrix(rep(0,p*q),nrow=p)
	zero2 <- matrix(rep(0,q*p),nrow=q)
	G3n <- rbind(  cbind(G31n,zero1)  , cbind(zero2,G32n)   )
	G5n <- rbind(  cbind(G1n,zero1)  , cbind(zero2,G2n)   )

	G1ninv <- solve(G1n)# indifinv(G1n,F)
	if(is.Fault(G1ninv)) G1ninv <- ginv(G1n)
	G2ninv<- solve(G2n)# indifinv(G2n,F)
	if(is.Fault(G2ninv)) G2ninv=ginv(G2n)
	G5ninv <- rbind(  cbind(G1ninv,zero1), cbind(zero2,G2ninv) )

	expr1 <- n * G5ninv
	expr2 <- (1/n) * G3n
	prd1 <- expr1 %*% expr2
	vcov <- prd1 %*% expr1
	vcov <- vcov / n
	parstdev = sqrt(diag(vcov[1:p,1:p]))				

	.expr4 <- sqrt((p+1)*qf(confidence,p+1,n-p-1))
	.expr5 <- parstdev * .expr4
	cilow <- unlist(object$parameters[names(object$form$par)]) - .expr5
	ciupp <- unlist(object$parameters[names(object$form$par)]) + .expr5
	result <- list(covmat=vcov[1:p,1:p],covtau=vcov[(p+1):(p+q),(p+1):(p+q)],parstdev=parstdev,CI=cbind(cilow,ciupp))
	return(result)
}


#**********************************************
