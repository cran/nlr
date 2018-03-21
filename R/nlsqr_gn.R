
#****************************************************************************************
#***                                                                                  ***

#***                           weight is vmat not the inverse                         ***
#***                                                                                  ***
#****************************************************************************************

nlsqr.gn <- function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.0010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),vm,rm=eiginv(t(chol(vm))) )
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda
	trace <- control$trace

	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
#	if(class(formula) != "nl.form") formula=nl.form(form=formula,par=start,independent=names(start))
	trial<-start
	datalist <- as.list(data)
	fact <- 10
	Fault2 <- Fault()
	p <- formula$p
  .parameters <- names(start)
  datalist[.parameters] <- start 												#*****datalist contains both parameter          **
	th <- start																		#*****vectors and data values                   **
	fmod <- eval(formula, datalist)												#***** First model computing                    **
	if(is.Fault(fmod)){
		qrfittresult <- fmod
		qrfittresult@FF<-"nsqr.gn"
		return(nl.fitt.gn(Fault=qrfittresult))
	}
	if(is.null(fmod)){																#** Can not compute the model, miss defined     **
		Fault2 <- Fault(FN=6,FF="nlsqr.gn")										#** model or bad choice for starting point      **
		return(nl.fitt.gn(Fault=Fault2))
	}
	n<-length(fmod$response)
	ndot <- n - p
	g <- attr(fmod$predictor, "gradient")
	h <- attr(fmod$predictor, "hessian")
	fmod$predictor<- transformNR(fmod$predictor,rm)
  fmod$response<- transformNR(fmod$response,rm)
	grd <- attr(fmod$predictor, "gradient")
	if (is.null(grd))
	{
		Fault2 <- Fault(FN=6,FF="nlsqr.gn")										#** The gradient attribute is not defined.      **
		return(nl.fitt.gn(Fault=Fault2))
	}
	else if (is.null(attr(fmod$predictor, "hessian")))
	{
		Fault2 <- Fault(FN=5,FF="nlsqr.gn")										#** The hessian attribute is not defined.       **
		return(nl.fitt.gn(Fault=Fault2))
	}
	mult <- sqrt(ndot / p)
	landa <- 1
	iteration <- 0
	iterhist <-NULL
	yresp <- as.numeric(fmod$response)		#y is response cane be any formula of y
	ybar <- mean(fmod$response)
															#*************************************************
	while (iteration <= maxiter)						#*********     Start Iteration    ****************
	{														#*********                        ****************
		iteration <- iteration + 1
		y <- as.numeric(fmod$response)					#y is response cane be any formula of y
		z0 <- y - as.numeric(fmod$predictor)
		dec <- qr(attr(fmod$predictor, "gradient"))
    if (dec$rank != p)
		{
			Fault2 <- Fault(FN=2,FL=F,FF="nlsqr.gn")	#** Singular Gradient matrix             ********		
			result <- nl.fitt.gn(
						parameters =   th,
						form =         formula,
						response =     transforminv(fmod$response,rm),
						predictor =    transforminv(fmod$predictor,rm), 
						history =      iterhist,
						method = 		   fittmethod(methodID=			11,
														methodBR=			14,
														detailBR=			"Modified Newton",
														subroutine=		"nlsqr.gn"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2,
						vm =           vm,
						rm =           rm,
						gresponse =    fmod$response,
						gpredictor =   fmod$predictor
					)				
			return(result)			
		}
		delta0 <- qr.coef(dec, z0)
		qtz0 <- qr.qty(dec, z0)
		ssq <- sum(z0 ^ 2)
		converge <- mult * sqrt(sum(qtz0[1:p] ^ 2)) / sqrt(sum(qtz0[n - (1:p)] ^ 2))
		if (converge < tolerance)
		{
			if(iteration==1) SumSquare <- ssq
			iterhist <- rbind(iterhist,c(iteration=iteration,objfnc=SumSquare, unlist(th), converge=converge))
			break
		}
		while (landa >= minlanda)   									 #========  Start iteraton for compute Landa  ======
		{                          										 #========                                    ======
			for (j in 1:p) trial[j] <- th[[j]] + landa * delta0[j]
      datalist[.parameters] <- trial
			fmod2 <- eval(formula, datalist)
			if(is.Fault(fmod2)){
				qrfittresult <- fmod2
				qrfittresult@FF<-"nsqr.gn"
				return(nl.fitt.gn(Fault=qrfittresult))
			}
			g <- attr(fmod2$predictor, "gradient")
			h <- attr(fmod2$predictor, "hessian")
      fmod2$predictor<- transformNR(fmod2$predictor,rm)
			fmod2$response<- transformNR(fmod2$response,rm)

      if(is.null(fmod2)){
				Fault2 <- Fault(FN=7,FF="nlsqr.gn")					#** Can not compute the model, missing mybe **				
				qrfittresult<-nl.fitt.gn(Fault=Fault2)
				return(qrfittresult)
			}
			z2 <- y - fmod2$predictor[1:n]
      ssq2 <- sum(z2 ^ 2)
			if (ssq2 <= ssq){
				break
			}
			landa <- landa / fact
		}                           									#==========    End Iteration Landa        ==========
		                          										#===================================================  
		if (landa >= minlanda)
		{
			SumSquare <- ssq2
			resid <- z2
		}
		else
		{
			resid <- z0
			SumSquare <- ssq
		}
		landa <- min(1, fact * landa)
		th <- trial                          						#**********   Last parameter values  ******
		datalist[.parameters] <- th
		fmod <- eval(formula, datalist)   							#**********   Last model computing   ******

    if(is.null(fmod)){
				Fault2 <- Fault(FN=7,FF="nlsqr.gn")				#** Can not compute the model, missing mybe **				
				qrfittresult<-nl.fitt.gn(Fault=Fault2)
				return(qrfittresult)
		}
		if (any(is.na(attr(fmod$predictor, "gradient"))))
		{
			Fault2 <- Fault(FN=6,FF="nlsqr.gn")					#** The gradient attribute is not defined.      **
			return(nl.fitt.gn(Fault=Fault2))
		}
		else if (any(is.na(attr(fmod$predictor, "hessian"))))
		{
			Fault2 <- Fault(FN=5,FF="nlsqr.gn")					#** The hessian attribute is not defined.       **
			return(nl.fitt.gn(Fault=Fault2))
		}
		g <- attr(fmod$predictor, "gradient")
		h <- attr(fmod$predictor, "hessian")
    fmod$predictor<- transformNR(fmod$predictor,rm)
		fmod$response<- transformNR(fmod$response,rm)		
    iterhist <- rbind(iterhist,c(iteration=iteration,objfnc=SumSquare, unlist(th), converge=converge))
	}      #*******************     Ent Iteration    *******************************
	       #************************************************************************
	if (iteration > maxiter) 
		Fault2 <- Fault(FN=1,FF="nlsqr.gn")				#** Can not compute the model, missing mybe **				
	sigma<-sqrt(iterhist[length(iterhist[,1]),"objfnc"]/(n-p))
	ssto <- sum( (as.numeric( fmod$predictor ) - ybar) ^ 2 )
	rho <- 1 - sum( ( yresp - as.numeric( fmod$predictor ) ) ^2 ) / ssto
	#curvatures<-curvature(attr(fmod$predictor,"gradient"),attr(fmod$predictor,"hessian"),sigma)
	curvatures<-NULL
  names(sigma) <- NULL
	th["sigma"] = sigma  # delete this latter, is done in scale
	tfm <- transforminv(fmod$predictor,rm)

  result <- nl.fitt.gn(parameters =  th,
            scale = sigma,
						correlation =  rho,
						form =         formula,
						response =     transforminv(fmod$response,rm),
						predictor =    tfm, 
						curvature =    curvatures,
						history =      iterhist,
						method =       fittmethod(methodID=			11,
											methodBR=			14,
											detailBR=			"Modified Newton",
											subroutine=		"nlsqr.gn"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2,
						vm =           vm,
						rm =           rm,
						gresponse =    fmod$response,
						gpredictor =   fmod$predictor
					)
  return(result)
	
#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#***                            End of File qrfit.                                    ***
#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
}
