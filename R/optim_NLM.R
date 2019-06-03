
#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |   Function optim.NLM, optimize function by levenberg marquardt.        | *
#* |                                                                        | *
#* |  Note: becarefull to using this function when there is not outlier, it | *
#* |    may not work witout outlier, in this case better to use nlmest      | *
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
optim.NLM<-function(objfnc,data,start=getInitial(objfnc,data),
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
		ilev <-0
		if(is.Fault(hshinv))  				#####################################  LM modified
		{
								###################################################################
								###      LM part
								###################################################################
												if(any(is.infinite(hsh)) || any(is.na(hsh))) return(Fault(FN=18))
												heig <-eigen(hsh,symmetric=T)         ## eig(f+land I) = (eig(F) +lnd)
												lambda <- abs(min(heig$values))*2  ## add the smallest minus eig to
                                     													## to make all possitive
												repeat{
													ilev<- ilev + 1
													I <- diag(dim(hsh) [2])
													zeta <- hsh + lambda * diag(abs(diag(hsh)))
													zetainv <- eiginv(zeta,stp=F,symmetric=T)
													if(is.Fault(zetainv))  return(nl.fitt.rob(Fault=zetainv))
														#zetainv <- svdinv(hsh)
													zqr <-qr(zeta)
										      delta2 <- grh  %*%  zetainv
													theta2 <- theta1 - delta2
													th2 <- as.list(theta2)
													names(theta2) <- names(theta1)
													names(th2) <- names(theta1)
													ht2 <- objfnc(data=data,start=th2,...)
													if(is.Faultwarn(ht2)){
													 return(ht2)}
													cnv <- sum((theta1-theta2)^2)
													diference <- abs(as.numeric(ht2$value)-as.numeric(ht$value))
													if(as.numeric(ht2$value) <= as.numeric(ht$value) || cnv < tolerance/10^10)
													{
														theta1 <- theta2
														ht <- ht2
														lambda <- lambda / 10
														break
														if(cnv < tolerance/10^10) break
													}
													else {
														lambda <- lambda * 10
														if( lambda > 1e12) lambda <- lambda / 20#return(Fault(FN=14,FF="optim.NLM"))
													}
												}
		}
							#####################     end of singular case                 ####################
							###################################################################################
		else {  			#####################   Possitive definit case  ###################################
						delta1 <- grh %*% hshinv #%*% t(grh)
						delta1 <- as.numeric(delta1)
						theta2 <- theta1 - landa*delta1
						th2 <- as.list(theta2)
						names(theta2) <- names(theta1)
						names(th2) <- names(theta1)
				        ht2 <- objfnc(data=data,start=th2,...)
						if(is.Faultwarn(ht2)) {
							 return(ht2)			
				            }			
						cnvg <- cnvgnew <- F
						diference <- abs(as.numeric(ht2$value)-as.numeric(ht$value))
#						cat("positive definit case",diference ,tolerance/1000,"\n")
#						cat("theta 1...theta2==",theta1,theta2,"\n","values 2..1=",as.numeric(ht2$value), as.numeric(ht$value))
#						cat("\n landa... delta",landa,delta1,landa*delta1,"\n \n")
						.temp1<-as.numeric(ht2$value)
						.temp2<-as.numeric(ht$value)
						if(is.missing(.temp1) ||	is.nan(.temp1) || is.inf(.temp1) ||
							is.missing(.temp1) ||	is.nan(.temp1) || is.inf(.temp1)) return(nl.fitt.rob(Fault=Fault(FN=18,FF="optim.NLM")))
						if(as.numeric(ht2$value) > as.numeric(ht$value)){
							ht.prev <- as.numeric(ht2$value)
														###############################
							while(landa >= minlanda){ 			####### iteration landa     ###
									landa <- landa / fact
									theta2 <- theta1 - landa * delta1
									th2 <- as.list(theta2)
									names(theta2) <- names(theta1)
									names(th2) <- names(theta1)
									ht2 <- objfnc(data=data,start=th2,...)
									if(is.Faultwarn(ht2)) return(ht2)
									ht.new <- as.numeric(ht2$value)
									diference <- abs(as.numeric(ht2$value)-as.numeric(ht$value))									
									if(as.numeric(ht2$value) <= as.numeric(ht$value))
										{
											theta1 <- theta2
											ht <- ht2
											cnvg <- T
											break
										}
									if(ht.new < ht.prev) {
										ht2.new <- ht2
										theta2.new <- theta2
										ht.prev <- ht.new
										cnvgnew <- T
									}
							}										#######  iteration landa    ###
																	###############################
							landa <- min(1,fact*landa)
							#landa <- 1
							if(! cnvg)
								if(cnvgnew){
									theta1 <- theta2.new
									ht <- ht2.new
								}
						}
						else{
							theta1 <- theta2
							ht <- ht2
							#landa <- landa * fact
							cnvg <- T
						}
		}			#######################     end positive                                      ############
					#######################  from above theta1 & ht must be returned back         ############
					##########################################################################################
		g2 <- as.matrix(ht$angvec)
		htg <- ht$angmat
		.expr1 <- t(htg) %*% (htg)											## H' H  p*p
		.expr2 <- eiginv(.expr1,symmetric=T,stp=F)						## (H' H)-1
		if(is.Fault(.expr2)){
			heig <-eigen(.expr1,symmetric=T)
			lambda <- abs(max(heig$values)) 
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
		if(is.missing(angle)||	is.nan(angle) || is.inf(angle)) return(Fault(FN=16,FF="optim.NLM"))
		th <-as.list(theta1)													##  without sigma
		names(th) <- names(theta1)
		iterhist <- rbind(iterhist,c(iteration = iterate,objfnc = as.numeric(ht$value),
											unlist(th),converge = angle	,ilev=ilev))
		if(angle < tolerance || diference < (tolerance/10)*abs(as.numeric(ht$value))) eol <- T
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





#******************************************************************************
#* +------------------------------------------------------------------------+ *
#* |                         End of 'optim.NLM'                             | *
#* |                                                                        | *
#* |                 Hossein Riazoshams, UPM, INSPEM                        | *
#* |                                                                        | *
#* |                          modified from 'nlmest'                        | *
#* |                                                                        | *
#* |                               11  Jan 2010                             | *
#* |                                                                        | *
#* +------------------------------------------------------------------------+ *
#******************************************************************************











