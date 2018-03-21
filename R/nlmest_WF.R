#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function 'nlmest', M-estimate of a nonlinear function.               | *
#* |       Using Wolf Conditions.                                           | *
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

nlmest.WF<-function(formula,data,start=getInitial(formula,data),robfunc,
	control=nlr.control(tolerance=0.0001, maxiter=25 * length(start),robscale=T),vm=NULL,
	rm=eiginv(t(chol(vm))),...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	loc.start <- start
	if(is.null(vm)) loc.objfnc <- robloss
	else  loc.objfnc <- robloss.gn
	if(trace) plot(data[[formula$independent]],data[[formula$dependent]])
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	if (is.null(loc.start$sigma)){
		dt1 <- c(data,loc.start)
		ht <- eval(formula,dt1)
		dv <- as.numeric(ht$predictor) - as.numeric(ht$response)
		if(! is.null(vm)) dv <- rm %*% dv
		loc.start[["sigma"]] <- mad(dv)
	}
#	print("loc.start[[sigma]] =====")
#	print(loc.start[["sigma"]] )
	th <- loc.start                                               ##  with sigma
	theta <- unlist(loc.start[names(loc.start)!="sigma"])         ##  without sigma
	theta1 <- theta                                               ##  without 
	datalist<-as.list(data)
	p <- length(theta)
	n <- length(data[[1]])
  	.parameters <- names(start)
	datalist[.parameters] <- loc.start[.parameters]      #*****datalist contains both parameter vectors and data values
	eol <- F
	iterate <- 0
	iterhist <- NULL
	lambda <-1
	if(is.null(vm)) 	ht <- robloss(formula,data,th,robfunc,control=control,...)		## th: with sigma
	else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
	tmp <- as.numeric(ht$ri)
	if(control$robscale) sigma<-mscale(tmp) # <- nl.mscale(tmp,robfunc,...)
	else sigma  <- sum(abs(tmp)) / n
	if(is.Fault(sigma)) return(sigma)
	if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
	if(trace) lines(data[[formula$independent]],as.numeric(ht$fmod$predictor))
	yresp <- as.numeric(ht$fmod$response)
	ybar <- mean(yresp)
	landa <- 1
	fact <- 2
						      #///*************************************************
	while (!eol)		#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
#   cat("===========================================iterate ================================= ",iterate,"\n")
		grh <- attr(ht$htheta,"gradient")									## g(theta)  grad
		hsh <- attr(ht$htheta,"hessian")									## H(theta)  hess
		hshinv <- eiginv(hsh,symmetric=T,stp=F)
		ilev <-0
		if(is.Fault(hshinv))  				#####################################  LM modified
		{
		###################################################################
		###      Negativ Def part
		###################################################################
				if(any(is.infinite(hsh)) || any(is.na(hsh))) return(Fault(FN=18))
				heig <-eigen(hsh,symmetric=T)         ## eig(f+land I) = (eig(F) +lnd)
				lambda <- 1  ## add the smallest minus eig to
            					## to make all possitive
				damping <-  diag(dim(hsh) [2])
				eigvh <-heig$values
				eigsum <- sum(abs(eigvh))
				eignorm <- eigvh/eigsum
				dltc=1e-8
				smalleig <- eigvh <dltc
				eigvh[smalleig]<-dltc-eigvh[smalleig]
				eigvh <- diag(eigvh)
				hmodif <- hsh + eigvh
				.temp2<-(heig$vectors) %*% eigvh
				zeta <- .temp2 %*% t(heig$vectors)
				hshinv<- eiginv(zeta,stp=F,symmetric=T)
				if(is.Fault(hshinv)) return(nl.fitt.rob(Fault=hshinv))
				######### cholesky with addeed multiple of the identity
#				beta <- 1e-3
#				if(min(diag(hsh)>0))  t0 <- beta
#				else t0 <- -min(diag(hsh)) + beta
#				repeat{
#					Imat <- diag(dim(hsh)[2])
#					zeta <- hsh + t0 * Imat
#					chsh <- chol(zeta)
#					if(attr(chsh,"rank")< p ) t0 <-max( 2*t0,beta)
#					else break
#				}
#				hshinv<- eiginv(zeta,stp=F,symmetric=T)
#				if(is.Fault(hshinv)) return(nl.fitt.rob(Fault=hshinv))
#			print("singularity problem solved")
		}
		###################################################################
		delta1 <- hshinv %*% grh
		pk <- -delta1
		landamax <- 1
		landa2 <- 1
		landa1 <- 0
		c1<-1.0e-8
		c2<-0.9
																#** NOTE: ht is old f() at xk **********************************
		phi1 <- phi0 <-as.numeric(ht$htheta)           #** phi1 change inside itteration but phi0 is fixed anywhere   **
		phid1 <- phid0 <- as.numeric(grh %c% (pk))  #** phid1 change inside itteration but phid0 is fixed anywhere **	
		ilev <- 1
#		if(trace) if(formula$p ==2){
		####################### phi(alpha) graph ##############
#			x=seq(-10,100,length.out=50)
#			y=NULL
#			z=NULL
#			derv=NULL
#			for(i in 1:length(x)){
#				t= theta1 - x[i]*delta1
#				loc.arg=as.list(t)
#				loc.arg["sigma"]<-sigma
#				if(is.null(vm)) t2 <- robloss(formula,data,loc.arg,robfunc,control=control,...)
#				else t2 <- robloss.gn(formula,data,loc.arg,robfunc,rm,control=control,...)
#				y[i]=as.numeric(t2$htheta)
#				z[i]=(phi0+x[i]*phid0)
#				derv[i]=as.numeric(attr(t2$htheta,"gradient") %c% (pk))     ## a%c%a is t(a) %*% a
#			}
#			plot(derv,type="l")
#			plot(x,y,type="l")
#			lines(x,z)
#			points(landa1,phi1)
#		}
		############################################################################
		########### now iterate, 
		########### 
		repeat{
#	      	print(theta1)
#     		 print(delta1)
			theta2 <- theta1 - landa2*delta1
#			cat("repppppp-------------------------------------------------------------------------------------------------------------------------  ilev=",ilev,"\n")
			loc.arg=as.list(theta2)
			names(loc.arg) <- names(formula$par)      
			loc.arg["sigma"]<-sigma			
			if(is.null(vm)) ht2 <- robloss(formula,data,loc.arg,robfunc,control=control,...)
			else ht2 <- robloss.gn(formula,data,loc.arg,robfunc,rm,control=control,...)
			if(is.Fault(ht2))  return(nl.fitt.rob(Fault=Fault(FN=18,FF="nlmest.WF")))
			.temp1<-as.numeric(ht2$htheta)
			diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/abs(as.numeric(ht$htheta)+.005)   ### can be used after iteration
			phi2 <- as.numeric(ht2$htheta)

			if(trace) points(landa2,phi2)
			ht2g <-attr(ht2$htheta,"gradient")
			phid2 <- as.numeric(ht2g %c% (-delta1))
#			cat("increament = ",-landa2*delta1,"\n")
#			cat(" ht2$value=",ht2$htheta,"ht$htheta= ",ht$htheta,"\n phi2= ",phi2,"(phi0+c1*landa2*phid0)",(phi0+c1*landa2*phid0),"\n (phid2) = ",(phid2),"-c2*phid0",-c2*phid0,"phid0 =",phid0, (phi2 > (phi0+c1*landa2*phid0)),"\n")			
#			cat("landa 2=", landa2,"theta2 =",theta2,"theta1 = ",theta1,"delta1 = ",delta1,"\n")
			if(  (phi2 > (phi0+c1*landa2*phid0)) || (phi2>phi1 & ilev>1)   ){     #**** zoom(alpha i-1 and alpha i) and stop
				zj=zoom2(landa1,landa2,phi1,phid2,phid1,phid2,ht,phi0,phid0,theta1=theta1,delta1=delta1,
				  sigma=sigma,objfnc=loc.objfnc,
					data=data,start=start,robfunc=robfunc,rmat=rm,formula=formula,control=control,...)
				if(is.Fault(zj)) return(zj)
				diference <- abs(as.numeric(zj$htj$htheta)-as.numeric(ht$htheta))/(as.numeric(ht$htheta)+.005)
				theta1 <- zj$thetaj
				ht<-zj$htj
				landa2 <- zj$landaj
				break
			}
			if(abs(phid2) <= -c2*phid0){
				diference <- abs(as.numeric(ht2$htheta)-as.numeric(ht$htheta))/(as.numeric(ht$htheta)+.005)
				theta1 <- theta2
				ht <- ht2
				break
			}
			if(phid2 >= 0) {  #**** zoom(alpha i and alpha i-1) and stop
				zj=zoom2(landa2,landa1,phi2,phid1,phid2,phid1,ht,phi0,phid0,theta1=theta1,delta1,sigma=sigma,objfnc=loc.objfnc,
					data=data,start=start,robfunc=robfunc,rmat=rm,formula=formula,control=control,...)
				if(is.Fault(zj)) return(zj)
				diference <- abs(as.numeric(zj$htj$htheta)-as.numeric(ht$htheta))/(as.numeric(ht$htheta)+.005)
				theta1 <- zj$thetaj
				ht<-zj$htj
				break
			}
			#*********************************
			.temp<-landa2
			landa2 <- landa2 + (landa1-landa2)*0.5   #runif(1,0,1)
#			zj=zoom2(landa1,landa2,phi1,phid2,phid1,phid2,ht,phi0,phid0,theta1=theta1,delta1,maxiter=maxiter,
#			sigma=sigma,trace=trace,objfnc=loc.objfnc,data=data,start=start,tolerance=tolerance,robfunc=robfunc,rmat=rm,control=control,...)
			landa1<-.temp
			ht<- ht2
			phi1<-phi2
			phid1<-phid2
			theta1<-theta2
			if(abs(landa2-landa1)<tolerance/1000) {
				break
			}
			ilev <- ilev+1
			if(ilev>maxiter) break
		}###############################
		##########################################################################################
		g2 <- attr(ht$rho,"gradient")									## V(heta) = [rho.(r1/sg)   .....  rho.(rn/sg)]T
		tmp <- as.numeric(ht$ri)
		if (control$robscale) sigma <-mscale(tmp)   # nl.mscale(tmp,robfunc,...)
		else sigma <- sum(abs(tmp)) / n
#		cat("\n sigma = ",sigma,"\n")
		if(is.Fault(sigma)) return(sigma)			
		.expr1 <- attr(ht$ri,"gradient") %c% attr(ht$ri,"gradient")	## J' J  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)
		if(is.Fault(.expr2)){
			heig <-eigen(.expr1,symmetric=T)
			lambda <- abs(min(heig$values)) * 2
			I <- diag(dim(.expr1) [2])
			.expr1 <- .expr1 + lambda * I
			.expr2 <- eiginv(.expr1,symmetric=T,stp=F)			
		}
		if(is.Fault(.expr2)) return(nl.fitt.rob(Fault=.expr2))
		.expr3 <- attr(ht$ri,"gradient") %*% .expr2						## J (J' J)^-1		n*p
		.expr4 <- .expr3 %*% t(attr(ht$ri,"gradient"))					## J (J' J)^-1 J'  n*n
		.expr5 <- .expr4 %*% g2												## VT = J (J' J)^-1 J' V 	n*1
		angle <- g2 %c% .expr5												## V' * VT		1*1
		angle <- angle / sqrt( sum(g2^2) * sum(.expr5^2) ) 
#		cat("angle   =",angle,"diference=",diference,"\n")
		th <-as.list(theta1)												##  without sigma
		names(th) <- names(formula$par)
		th["sigma"]  <- sigma												##  renew with sigma
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$htheta),
			unlist(th),converge = angle	,ilev=ilev))
		if(is.missing(angle)||	is.nan(angle) || is.inf(angle)) return(nl.fitt.rob(Fault=Fault(FN=16,FF="nlmest.WF")))
		if(angle < tolerance || diference < tolerance/1000) eol <- T 
		else 
			if(iterate > maxiter) {
				eol <- T
				Fault2 <- Fault(FN=1,FF = "nlmest")
			}
			else {
				if(is.null(vm)) ht <- robloss(formula,data,th,robfunc,control=control,...)  ## h(theta) = sum(  rho(ri/sigma)  )
				else ht <- robloss.gn(formula,data,th,robfunc,rm,control=control,...)
				if(trace) lines(data[[formula$independent]],ht$fmod$predictor)
#				cat("222222 sigma= ",sigma,"\n")
				if(ht$Fault$FL) return(nl.fitt.rob(Fault=ht$Fault))
				tmp <- as.numeric(ht$ri)
				if (control$robscale) sigma <-mscale(tmp)   # nl.mscale(tmp,robfunc,...)
				else sigma <- sum(abs(tmp)) / n
				th["sigma"]  <- sigma												##  renew with sigma				
			}
			#\\\*********                        **********************************************************************************
	}		#\\\*********     End Iteration      **********************************************************************************
			#\\\*******************************************************************************************************************
	#//////////*******************************************************************************************************************
	#//////////*******************************************************************************************************************
	#//////////********     preparing output   ***********************************************************************************
	#//////////********                        ***********************************************************************************
	if(! is.null(vm)){
		rinv <- eiginv(rm)
		#yresp <- rinv %*% as.numeric(ht$fmod$response)
		#ypred <- rinv %*% as.numeric(ht$fmod$predictor)
		yresp <- as.numeric(ht$fmod$response)
		ypred <- as.numeric(ht$fmod$predictor)		
	}
	else{
		yresp <- as.numeric(ht$fmod$response)
		ypred <- as.numeric(ht$fmod$predictor)
	}
	ybar <- mean(yresp)	
	nlrho <- 1 - sum( ( yresp - ypred )  ^ 2 ) / 
					sum( (ypred - ybar ) ^ 2 )
	curv1 <- curvature(gradient = attr(ht$fmod$predictor,"gradient"),
						hessian = attr(ht$fmod$predictor,"hessian"),
						sigma = th[["sigma"]])

  if(is.null(vm))
		result <- nl.fitt.rob(
            parameters =  th,
            scale = sigma,
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
						predictor =    ht$fmod$predictor,
						curvature =    curv1,
						history =      iterhist,
						method =       fittmethod(methodID=			7,
				                      methodBR=			switch(control$robscale,10,11),       ### rob or nopnrob variance
														  detailBR=			"Lev-Marq",
														  subroutine=		"nlmest.WF"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =      robfunc
					)
	else
		result <- nl.fitt.rgn(
            parameters =  th,
            scale = sigma,
						correlation =  nlrho,
						form =         formula,
						response =     ht$fmod$response,
            predictor =    ht$fmod$predictor, 
						curvature =    curv1,
						history =      iterhist,
						method = 		    fittmethod(methodID= 8,
													  	methodBR=			switch(control$robscale,10,11),       ### rob or nopnrob variance
													  	detailBR=			"Lev-Marq",
													  	subroutine=		"nlmest.WF"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2,
            others =       NULL,
						htheta =       ht$htheta,
						rho =          ht$rho,
						ri =           ht$ri,
						curvrob =      NULL,
						robform =       robfunc,
						vm =           vm,
						rm=            rm,
						gresponse =    transformNR(ht$fmod$response,rm),
						gpredictor =   transformNR(ht$fmod$predictor,rm)
					)
	return(result)
}
   #****************************************************
   #**      Cubic interpolation                       ** 
   #**      Defined in optim.WF also                  **
   #****************************************************

   #****************************************************
   #**     end  Cubic interpolation                   ** 
   #****************************************************
   #****************************************************
   #**      zoom function                             ** 
   #**                                                **
   #****************************************************
zoom2<-function(a1,a2,p1,p2,pd1,pd2,ht,phi0,phid0,theta1,delta1,sigma,objfnc,data,start,control,...){
		c1<-1.0e-4
		c2<-0.9
		jlev<-1
					repeat{
						d1<-pd1+pd2-3*( (p1-p2)/(a1-a2) )       ### this cubic interpolation doesnt work well
						d2<- sign(a2-a1)* sqrt( d1^2 - pd1*pd2)
						aj<- a2  -  (a2-a1) * (  (pd2+d2-d1)/(pd2-pd1+2*d2)  )
						#aj<-1/2*(a1+a2)
						thetaj <- theta1 - aj*delta1
						loc.arg<-as.list(thetaj)
            names(loc.arg) <- names(start[names(start) != "sigma"])
						loc.arg["sigma"] <- sigma
#						cat("increament = ",- aj*delta1,"\n")
#						cat("a1= ",a1,"a2= ",a2,"aj = ",aj,"\n theta j= ",thetaj,"theta 1= ",theta1,"\n")
#						cat("pd1= ",pd1,"pd2= ",pd2,"\n")
						htj <- objfnc(data=data,start=loc.arg,...)
						if(is.Fault(htj)){
							aj<- a2- (pd2*(a2-a1))/(pd2-pd1)         ### this is parabolic search see Gaviano 1987
							thetaj <- theta1 - aj*delta1
							thj<-as.list(thetaj)
							thj$sigma<-loc.arg$sigma
							names(thj) <- names(loc.arg)
							htj <- objfnc(data=data,start=thj,...)							
							if(is.Fault(htj)) return(htj)
						}
						pj <- as.numeric(htj$htheta)                 ## phi(aj)
						grhj<- attr(htj$htheta,"gradient")           ## f'(aj)
						pdj<- as.numeric(grhj %c% (-delta1))  		 ## phi'(aj)
						if(abs(aj-a2)<control$tolerance/1000) return(list(thetaj=thetaj,landaj=aj,htj=htj,pj=pj,pdj=pdj)) 
#						cat("phij= ", pj,"(phi0 + c1*aj*phid0)",(phi0 + c1*aj*phid0),"\n","(pdj)= ",(pdj),"(-c2*phid0)",(-c2*phid0),"pdj*(a2-a1)=",pdj*(a2-a1),"\n \n")
						if( (pj > (phi0 + c1*aj*phid0)) || (pj >= p1)){
							a2 <- aj
							p2 <- pj
							pd2 <- pdj
						}
						else{
							if(abs(pdj) <= (-c2*phid0)){
								return(list(thetaj=thetaj,landaj=aj,htj=htj,pj=pj,pdj=pdj)) 
							}
							if(pdj*(a2-a1) >=0) {
								a2<-a1
								p2<-p1
								pd2<-pd1
							}
							a1<-aj
							p1<-pj
							pd1<-pdj
						}
						jlev<-jlev+1
						if(jlev>control$maxiter) return(list(thetaj=thetaj,landaj=aj,htj=htj,pj=pj,pdj=pdj)) 
					}
}
   #****************************************************
   #**     end  Cubic interpolation                   ** 
   #****************************************************

#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of optim.WF'                               | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, Stockholm                          | *
#* |                                                                        | *
#* |                            13  Nov 2012                                | *
#* |                                                                        | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************









