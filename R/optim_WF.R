#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |  Function optim.WF, optimize function by Leaneat search Newton.        | *
#* |  minimizing alpha by quadratic interpolation, and check WOLF conditions| *
#* |                                                                        | *
#* |    argumnts:                                                           | *
#* |      objfnc: any objective function for minimizing, it must contains:  | *
#* |              formula, data and start, extra will be defined in (...)   | *
#* |              the output of objfnc must contains:                       | *
#* |                $value(attr,gradient,hessian), $angmat,$angvec          | *
#* |      data:    data, contains dependents and independents,              | *
#* |               data.frame, list or named matrix it can be.              | *
#* |      theta:   starting values, it must contains tau elements           | *
#* |                                                                        | *
#* |      ...:     can be entries objfnc                                    | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************
optim.WF<-function(objfnc,data,start=getInitial(objfnc,data),
	control=nlr.control(tolerance=0.001, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),...)
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	loc.start <- start	
	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
	th <- loc.start
	theta <- unlist(loc.start)
	theta1 <- theta
	p <- length(theta)
	n <- length(data[[1]])
	eol <- F
	iterate <- 0
	iterhist <- NULL
	.datalist <- NULL
	.datalist[names(loc.start)] <- loc.start
	.datalist[names(data)]<-data
	ht <- objfnc(data=data,start=th,...)
	if(is.Faultwarn(ht)) return(ht)
	lambda <- 1
	landa <-1
	fact <- 2
						#///*************************************************
	while (!eol)		#///*********     Start Iteration    ****************
	{
		iterate <- iterate + 1
		grh <- attr(ht$value,"gradient")                          ## g(theta)  grad
		hsh <- attr(ht$value,"hessian")                           ## H(theta)  hess
		hshinv <- eiginv(hsh,stp=F,symmetric=T)
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
		}
		###################################################################
		delta1 <- hshinv %*% grh
		landamax <- 1
		landa2 <- 1-1e-8
		landa1 <- 0
		c1<-1.0e-8
		c2<-0.9
															#** NOTE: ht is old f() at xk **********************************
		phi1 <- phi0 <-as.numeric(ht$value)           #** phi1 change inside itteration but phi0 is fixed anywhere   **
		phid1 <- phid0 <- as.numeric(t(grh) %*% (-delta1))  #** phid1 change inside itteration but phid0 is fixed anywhere **	
		ilev <- 1
		if(trace){
		####################### phi(alpha) graph ##############
			x=seq(-10,100,length.out=50)
			y=NULL
			z=NULL
			derv=NULL
			for(i in 1:length(x)){
				t= theta1 - x[i]*delta1
				t=as.list(t)
				t2 <- objfnc(data=data,start=t,...)
				y[i]=as.numeric(t2$value)
				z[i]=(phi0+x[i]*phid0)
				derv[i]=as.numeric(t(attr(t2$value,"gradient")) %*% (-delta1))
			}
			plot(derv,type="l")
			plot(x,y,type="l")
			lines(x,z)
			points(landa1,phi1)
		}
		############################################################################
		########### now iterate, 
		########### 
#		print("goooooooooooooooooooooooooooooooooooooooooooooooooooo")
		repeat{
			theta2 <- theta1 - landa2*delta1
			th2 <- as.list(theta2)
			ht2 <- objfnc(data=data,start=th2,...)
			.temp1<-as.numeric(ht2$value)
			if(is.Fault(ht2)) return(Fault=Fault(FN=18,FF="optim.WF"))
			diference <- abs(as.numeric(ht2$value)-as.numeric(ht$value))/(as.numeric(ht$value)+.005)	 ### can be used after iteration			
			phi2 <- as.numeric(ht2$value)

			if(trace) points(landa2,phi2)
			ht2g <-attr(ht2$value,"gradient")
#			cat("reppppppppppppppppppppppppppppppp  ilev=",ilev,"\n")
			phid2 <- as.numeric(t(ht2g) %*% (-delta1))
#			cat(" ht2$value=",ht2$value,"ht$value= ",ht$value,"theta 2 =",theta2,"\n","phi2= ",phi2,"phi0=",phi0,"(phi0+c1*landa2*phid0)",(phi0+c1*landa2*phid0),"(phid2) = ",(phid2),"-c2*phid0",-c2*phid0,"phid0 =",phid0, (phi2 > (phi0+c1*landa2*phid0)),"\n")			
#			cat("landa 2=", landa2,"theta2 =",theta2,"theta1 = ",theta1,"delta1 = ",delta1,"\n")

			if(  (phi2 > (phi0+c1*landa2*phid0)) || (phi2>phi1 & ilev>1)   ){     #**** zoom(alpha i-1 and alpha i) and stop
				zj=zoom(landa1,landa2,phi1,phid2,phid1,phid2,ht,phi0,phid0,theta1=theta1,delta1,maxiter=maxiter,objfnc=objfnc,data,start,trace=trace,tolerance=tolerance,...)
				if(is.Fault(zj)) return(zj)
				diference <- abs(as.numeric(zj$htj$value)-as.numeric(ht$value))/(as.numeric(ht$value)+.005)
				theta1 <- zj$thetaj
				ht<-zj$htj
				landa2 <- zj$landaj
				break
			}
			if(abs(phid2) <= -c2*phid0){
				diference <- abs(as.numeric(ht2$value)-as.numeric(ht$value))/(as.numeric(ht$value)+.005)
				theta1 <- theta2
				ht <- ht2
				break
			}
			if(phid2 >= 0) {  #**** zoom(alpha i and alpha i-1) and stop
				zj=zoom(landa2,landa1,phi2,phid1,phid2,phid1,ht,phi0,phid0,theta1=theta1,delta1,maxiter=maxiter,objfnc=objfnc,data,start,trace=trace,tolerance=tolerance,...)
				if(is.Fault(zj)) return(zj)
				diference <- abs(as.numeric(zj$htj$value)-as.numeric(ht$value))/(as.numeric(ht$value)+.005)
				theta1 <- zj$thetaj
				ht<-zj$htj
				break
			}
			#*********************************
			.temp<-landa2
			landa2 <- landa2 + (landa1-landa2)*.5 #runif(1,0,1)

			#zj=zoom(landa1,landa2,phi1,phi2,phid1,phid2,ht,phi0,phid0,theta1=theta1,
			#			delta1,maxiter=maxiter,objfnc=objfnc,data,start,trace=T,...)
			landa1<-.temp
			ht<- ht2
			phi1<-phi2
			phid1<-phid2
			theta1<-theta2
			if(abs(landa2-landa1)<tolerance/1000) break
			#landa2 <- zj$landaj
			#phi2<- zj$pj
			#phid2<- zj$pdj
			ilev <- ilev+1
			if(ilev>maxiter) break
		}###############################
		##########################################################################################
		g2 <- as.matrix(ht$angvec)
		htg <- ht$angmat
		.expr1 <- t(htg) %*% (htg)											## H' H  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)						## (H' H)-1
		if(is.Fault(.expr2)){
			heig <-eigen(.expr1,symmetric=T)
			lambda <- abs(max(heig$value)) 
			I <- diag(dim(.expr1) [2])
			.expr1 <- .expr1 + lambda * I
			.expr2 <- eiginv(.expr1,symmetric=T,stp=F)			
		}		
		#.expr2 <- zetainv
		if(is.Fault(.expr2)){
		 return(.expr2) 								## Note: eiginv only shoud return back Fault
		}
		.expr3 <- htg %*% .expr2                                      ## H (H' H)^-1		n*p
		.expr4 <- .expr3 %*% t(htg)											## H (H' H)^-1 H'  n*n
		.expr5 <- .expr4 %*% g2												## VT = H (H' H)^-1 H' V 	n*1
		angle <- t(g2) %*% (.expr5)											## V' * VT		1*1
		angle <- angle / sqrt( sum(g2^2) * sum(.expr5^2) )				## cos(V,VT)
		th <-as.list(theta1)													##  without sigma
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$value),
											unlist(th),converge = angle	,ilev=ilev))
		if(is.missing(angle)||	is.nan(angle) || is.inf(angle)) return(Fault(FN=16,FF="optim.WF"))
		if(angle < tolerance || diference < (tolerance/1000)) eol <- T
		else if(iterate > maxiter) {
			eol <- T
			Fault2 <- Fault(FN=1,FF = "optim.NLM")
		}
		else {
			ht <- objfnc(data=data,start=th,...)
			if(is.Faultwarn(ht)){
			 return(ht)
			}			
		}
	}		#\\\*********     End Iteration      ****************
			#\\\*************************************************
	result=list(parameters = th, objfnc=ht, history=iterhist)
	return(result)
}



	#****************************************************
	#*****   output:                                ***** 
	#*****  Structure of minus of pseudo likelihood ***** 
   #**                                                **
   #**    (very important) the output is different    **
   #**               from standard                    **
   #**  inheritance: nl.fitt:                         **
   #**      parameters = tau                          **
   #****************************************************
   #****************************************************
   #**      Cubic interpolation                       ** 
   #**                                                **
   #****************************************************
CubInterp<-function(a1,a2,p1,p2,pd1,pd2){
	d1<-pd1+pd2-3*( (p1-p2)/(a1-a2) )
	d2<- sign(a2-a1)* sqrt( d1^2 - pd1*pd2)
	result<- a2  -  (a2-a1) * (  (pd2+d2-d1)/(pd2-pd1+2*d2)  )
	result
}
   #****************************************************
   #**     end  Cubic interpolation                   ** 
   #****************************************************
   #****************************************************
   #**      zoom function                             ** 
   #**                                                **
   #****************************************************
zoom<-function(a1,a2,p1,p2,pd1,pd2,ht,phi0,phid0,theta1,delta1,maxiter=25 * length(theta1),objfnc,data,start,trace=F,tolerance,...){
		c1<-1.0e-4
		c2<-0.9
		jlev<-1
					repeat{
						d1<-pd1+pd2-3*( (p1-p2)/(a1-a2) )       ### this cubic interpolation doesnt work well
						d2<- sign(a2-a1)* sqrt( d1^2 - pd1*pd2)
						aj<- a2  -  (a2-a1) * (  (pd2+d2-d1)/(pd2-pd1+2*d2)  )
						#aj<-1/2*(a1+a2)
						thetaj <- theta1 - aj*delta1
						thj<-as.list(thetaj)
#						cat("a1= ",a1,"a2= ",a2,"aj = ",aj,"theta j= ",thetaj,"theta 1= ",theta1,"\n")
#						cat("pd1= ",pd1,"pd2= ",pd2,"\n")
						htj <- objfnc(data=data,start=thj,...)
						if(is.Fault(htj)){
							aj<- a2- (pd2*(a2-a1))/(pd2-pd1)         ### this is parabolic search see Gaviano 1987
							thetaj <- theta1 - aj*delta1
							thj<-as.list(thetaj)
							htj <- objfnc(data=data,start=thj,...)							
							if(is.Fault(htj)) return(htj)
						}
						pj <- as.numeric(htj$value)                 ## phi(aj)
						if(trace) points(aj,pj)
						grhj<- attr(htj$value,"gradient")            ## f'(aj)
						pdj<- as.numeric(t(grhj) %*% (-delta1))   ## phi'(aj)
						if(abs(aj-a2)<tolerance/1000) return(list(thetaj=thetaj,landaj=aj,htj=htj,pj=pj,pdj=pdj)) 
#						cat("phij= ", pj,"(phi0 + c1*aj*phid0)",(phi0 + c1*aj*phid0),"\n","abs(pdj)= ",abs(pdj),"(-c2*phid0)",(-c2*phid0),"pdj*(a2-a1)=",pdj*(a2-a1),"\n")
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
						if(jlev>maxiter) return(list(thetaj=thetaj,landaj=aj,htj=htj,pj=pj,pdj=pdj)) 
					}
}
   #****************************************************
   #**     end  Cubic interpolation                   ** 
   #****************************************************
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of 'optim.NLM'                             | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, Stockholm                          | *
#* |                                                                        | *
#* |                            07  Nov 2012                                | *
#* |                                                                        | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************











