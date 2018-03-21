#******************************************************************************************************
#**		nlsqr: Gss Newton linear approximation, using QR decomposition                               **
#**		Fault variab:                                                                                **
#**		FN  =  1:  i.e. Maximum number of itteration Exceeded, The computedarameters are unreliable  **
#**			FL  =  F, Warning not Error.                                                             **
#**		FN  =  2:  i.e. Singular Gradient maix, the parameter can not be computed                    **
#**			FL  =  T, Error                                                                          **
#**		FN  =  3:  i.e. Landa redces bellow minimum, Convergce may not acheived, Weakconvergence.	 **
#**			FL  =  F,May be weak convergence Not Error.                                              **
#**		FN  =  4:  The gradient attribute is not defined.                                            **
#**			FL  =  T, An error.                                                                      **
#**		FN  =  5:  The hessian attribute is not defined.                                             **
#**			FL  =  T, An error.                                                                      **
#**		FN  =  6:  The formula is not correctly defined.                                            **
#**			FL  =  T, An error.                                                                      **
#**		FN  =  7:  Can not compute the model, miss defined model or bad choice for starting point.   **
#**			FL  =  T, An error.                                                                      **
#******************************************************************************************************

nlsqr <- function(formula, data, start=getInitial(formula,data), 
		control=nlr.control(tolerance=0.00010, minlanda=1 / 2 ^ 10, maxiter=25 * length(start)))
{

	tolerance <- control$tolerance
	maxiter <- control$maxiter
	minlanda <- control$minlanda

	Fault2 <- Fault()		###### no error, automatically will be created by the Fault object
#	if(class(formula) != "nl.form") formula=nl.form(form=formula,par=start,independent=names(start))
	trial<-start
	datalist <- as.list(data)
	fact <- 10
	Fault2 <- Fault()
  .parameters <- names(start)
	p <- formula$p
	datalist[.parameters] <- start 												#*****datalist contains both parameter          **
	th <- start																		#*****vectors and data values                   **
	fmod <- eval(formula, datalist)												#***** First model computing                    **

	if(is.Fault(fmod)){																#** Can not compute the model, miss defined     **
		Fault2 <- Fault(FN=6)														#** model or bad choice for starting point      **
		return(nl.fitt(Fault=Fault2))
	}
	n<-length(fmod$response)
	ndot <- n - p

	grd <- attr(fmod$predictor, "gradient")
	if (is.null(grd))
	{
		Fault2 <- Fault(FN=4)														#** The gradient attribute is not defined.      **
		return(nl.fitt(Fault=Fault2))
	}
	else if(any(is.nan(grd))) return(Fault(FN=4))
	else if (is.null(attr(fmod$predictor, "hessian")))
	{
		Fault2 <- Fault(FN=5)														#** The hessian attribute is not defined.       **
		return(nl.fitt(Fault=Fault2))
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
			Fault2 <- Fault(FN=2,FL=F)					#** Singular Gradient matrix             ********		
			result <- nl.fitt(
						parameters =   th,
						form =         formula,
						response =     fmod$response,
						predictor =    fmod$predictor, 
						history =      iterhist,
						method = 		 fittmethod(methodID=			3,
														methodBR=			14,
														detailBR=			"Modified Newton",
														subroutine=		"nlsqr2"),
						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2
					)				
			return(result)			
			d <-  attr(fmod$predictor, "gradient") %c% z0
			vchol <- chol(attr(fmod$predictor, "gradient") %c% attr(fmod$predictor, "gradient"))
			zd <- eiginv(vchol) %*% d
			delta0 <- solve(t(vchol),zd)
		}
		else delta0 <- qr.coef(dec, z0)
		qtz0 <- qr.qty(dec, z0)
		ssq <- sum(z0 ^ 2)
		converge1 <- mult * sqrt(sum(qtz0[1:p] ^ 2)) 
    converge2 <- sqrt(sum(qtz0[n - (1:p)] ^ 2))
    converge3 <- converge2 *  tolerance
    converge <- converge1 / converge2
		#if (converge1 <= converge3)
		#converge <- mult * sqrt(sum(qtz0[1:p] ^ 2)) / sqrt(sum(qtz0[n - (1:p)] ^ 2))
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
			if(is.null(fmod2) || is.Fault(fmod2) ){
				Fault2 <- Fault(FN=7)									#** Can not compute the model, missing mybe **				
				qrfittresult<-nl.fitt(Fault=Fault2)
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
				Fault2 <- Fault(FN=7)								#** Can not compute the model, missing mybe **				
				qrfittresult<-nl.fitt(Fault=Fault2)
				return(qrfittresult)
		}
		if (any(is.na(attr(fmod$predictor, "gradient"))))
		{
			Fault2 <- Fault(FN=6)														#** The gradient attribute is not defined.      **
			return(nl.fitt(Fault=Fault2))
		}
		else if (any(is.na(attr(fmod$predictor, "hessian"))))
		{
			Fault2 <- Fault(FN=5)														#** The hessian attribute is not defined.       **
			return(nl.fitt(Fault=Fault2))
		}

		iterhist <- rbind(iterhist,c(iteration=iteration,objfnc=SumSquare, unlist(th), converge=converge))
	}      #*******************     Ent Iteration    *******************************
	       #************************************************************************
	if (iteration > maxiter) 
				Fault2 <- Fault(FN=1)									#** Can not compute the model, missing mybe **				
	sigma<-sqrt(iterhist[length(iterhist[,1]),"objfnc"]/(n-p))
	ssto <- sum( (as.numeric( fmod$predictor ) - ybar) ^ 2 )
	rho <- 1 - sum( ( yresp - as.numeric( fmod$predictor ) ) ^2 ) / ssto
	curvatures<-curvature(attr(fmod$predictor,"gradient"),attr(fmod$predictor,"hessian"),sigma)
	names(sigma) <- NULL
  th["sigma"] = sigma          # delete this latter, scale added
  th <- as.list(th)
  result <- nl.fitt(parameters =  th,
						scale = sigma,
            correlation =  rho,
						form =         formula,
						response =     fmod$response,
						predictor =    fmod$predictor, 
						curvature =    curvatures,
						history =      iterhist,
						method =       fittmethod(methodID=3,
                                      methodBR=14,
						                          detailBR="Modified Newton",
                                      subroutine="nlsqr2"),

						data =         as.list(data),
						sourcefnc =    match.call(),
						Fault =        Fault2
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
nlsqr2 <- nlsqr
