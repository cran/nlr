

#+#################################################################################+
#**       nlgn.fitt:  a generalized nonlinear fitted method. for storing fitted vlues  |
#**				Inherit from nl.fitt with all Methods.                        /
#**       generalized form is (y-f)' V (y-f), where v(ei) = sigma^2 * V, is variance |
#**                      covariance matrix of residuals.                             |
#**			Extra Slots                                                       /
#**             vmat:    matrix part of var cov matrix of residuals.               |
#**             rmat:    choleskey decomposition of V=U'U, R = (U')^-1, this       |
#**                          choleskey decomposition is used to transfer both side  \
#**                          of model by  R*y = R*f + e, which is const varian      /
#**                                                                  --------------|
#+--------------------------------------------------+
#|      class:  nl.fitt.rob                         |
#+--------------------------------------------------+

setClass(Class="nl.fitt.gn",representation(
									vm =          "matrix",
									rm =          "matrix",
									hetro =       "nl.fittorNULL",
									autcorr =     "listorNULL",      ## the list can be assigned 
																	## directly to a time series fitt
									autpar =      "listorNULL",      ##  looks like repiting again some
																	## part of above time series fitt, but ok 
									gresponse =   "vectororMatrix",
									gpredictor =  "vectororMatrix"
									),
				 contains = "nl.fitt",
         prototype=prototype(autcorr=NULL,autpar=NULL)
		)
###################################################
##   constructor, will call initialize method
##       By default the RY and RF will be computed
###################################################
nl.fitt.gn<-function(parameters,scale=NULL ,correlation=NULL,form ,response =NULL,predictor ,
						curvature =NULL,history =NULL,method=NULL,data,sourcefnc=NULL,Fault=Fault(),others=NULL,
						vm,rm,      #### compulsury
						hetro=NULL,
						autcorr = NULL,autpar=NULL,
						gresponse = transformNR(response,rm),
						gpredictor = transformNR(predictor,rm)
						){
	if(is.Fault(Fault)) return(new("nl.fitt.gn",Fault=Fault))

  result <- new("nl.fitt.gn",
						nl.fitt(
								parameters =   parameters ,
                scale =        scale,
								correlation=   correlation,
								form =         form, 
								response =     response , 
								predictor =    predictor , 
								curvature =    curvature ,
								history =      history , 
								method =       method,
								data =         data,
								sourcefnc =    sourcefnc,
								Fault =        Fault,
								others =       others),
							vm =           vm,
							rm =           rm,
							hetro =        hetro,
							autcorr =      autcorr,
							autpar =       autpar,
							gresponse =    gresponse,
							gpredictor =   gpredictor
				)

  return(result)
}

###################################################
##   Method residuals, nl.fitt.gn
###################################################
setMethod(f="residuals",signature="nl.fitt.gn",
	definition=function(object,data=NULL,...){
		if(is.null(data)){
			rs <- as.numeric(object@response)
			pr <- as.numeric(object@predictor)
			rsd <- (rs-pr)
			grsd <- object@rm %*% rsd
		}
		else{
			datalist <- as.list(data)
			datalist[names(object@form@par)] <- object@parameters[names(object@form@par)]
			fmod <- eval(object@form,datalist)
			rs <- as.numeric(fmod@response)
			pr <- as.numeric(fmod@predictor)
			rsd <- (rs-pr)
			grsd <- object@rm %*% rsd
		}
	result <- structure(.Data=rsd,"gresiduals"=grsd)
	return(result)
	}
)


#######################################################
##   Method recalc, nl.fitt.gn, recalculate object "object
##   for new values "colin".  added 14/06/2010
#######################################################
setMethod(f="recalc",signature="nl.fitt.gn",
	definition=function(object="nl.fitt",...){
		objcall <- object@sourcefnc
		if(is.null(objcall)) return(Fault(FN=13,FF="recalc nl.fitt.gn"))
		dtarg <- names(objcall)=="data"
		objcall[dtarg] <- call("objcall")
		parinit <- object@parameters
		pararg <- names(objcall)=="start"
		objcall[pararg] <- call("parinit")
		if(! is.null(object@rm)){
			parinit <- object@vm
			pararg <- names(objcall)=="vm"
			objcall[pararg] <- call("parinit")
		}
		renew <- eval(objcall )
		return(renew)	
	}
)

###################################################
##   Method parInfer, nl.fitt.rob, parameter Infrences
###################################################
setMethod("parInfer","nl.fitt.gn",
	 function(object,confidence = .95){
			.expr1 <- attr(object@gpredictor,"gradient")							##  RF
      .expr2 <- t(.expr1) %*% .expr1										##  (RF)' RF
			covmatrix <- object@scale ^ 2 * indifinv( .expr2 )	##  sgm2 * (RF' RF)^-1 
			v2inv <- diag(1/sqrt(diag(covmatrix)))
			corrmat <- v2inv %*% covmatrix %*% v2inv
			parstdev = sqrt(diag(covmatrix ))

			n <- nrow(.expr1)
			p <- object$form$p
			.expr4 <- sqrt((p+1)*qf(confidence,p+1,n-p-1))
			.expr5 <- parstdev * .expr4
			cilow <- unlist(object$parameters[names(object$form$par)]) - .expr5
			ciupp <- unlist(object$parameters[names(object$form$par)]) + .expr5

			result <- list(covmat=covmatrix ,corrmat = corrmat,parstdev=parstdev,CI=cbind(cilow,ciupp))

			return(result)
	}
)

#######################################################
##   Method PI, nl.fitt, prediction Interval
##     Added to this object att 27/11/2012
#######################################################
setMethod(f="predictionI",signature="nl.fitt.gn",
	definition=function(nlfited,confidence=.95,data=NULL){
		if(is.null(data)){
			yhat <- predict(nlfited)
			data <-nlfited@data[[nlfited@form$independent]]
			grdy <- attr(yhat,"gradient")
			yhat <- as.numeric(yhat)
			m <- nrow(grdy)
			varmod <-  predict(nlfited$hetro)
			gt <- as.numeric(varmod)
		}
		else{
			yhat <- predict(nlfited,newdata=data)
			grdy <- attr(yhat,"gradient")
			yhat <- as.numeric(yhat)
      if(is.null(nlfited$data[[nlfited$hetro$form$independent]])) data[[nlfited$hetro$form$independent]] <- yhat ### H(f,tau)
      else data[[nlfited$hetro$form$independent]] <- data[[nlfited$form$independent]]
			varmod <- predict(nlfited$hetro,newdata=data)
			gt <- as.numeric(varmod)
			m <- nrow(grdy)
		}
		sigma <- nlfited$scale
		gpg<-parInfer(nlfited)$covmat             ## sigma embded here
		gpgi <- indifinv(gpg)

		t6<-rep(0,times=m)
		for(i in 1:m){
			t5<-grdy[i,]%*%gpg			#*** 1*p
			t6[i]<-t5%*%(grdy[i,])       #*** 1*1
		}
		rf <- qt(1-confidence,m-nlfited@form@p)
		if(! is.null(nlfited$hetro)){ # correct in new version for autocorrelated
		  vp <- gt+t6                              ## sigma is embded inside both of variance parts
		}
		else{
		  vp <- t6                              ## sigma is embded inside both of variance parts
		}
		svp <- sqrt(vp)
		up<-yhat+rf*svp
		lw<-yhat-rf*svp
		return(list(lowerbound=lw,upperbound=up))
	}
)
#'***********************************************************
#'*   Method atypicals, nl.fitt.gn, parameter Infrences
#'     Added to this object att 19/6/2014
#'*     corrected at 8/2/2017
#'*  According to riaz2009 et al only studres and elliptnorm
#'*     can detect outliers. 
#'*     Other measures shown here for feature development.
#'***********************************************************

setMethod(f="atypicals",signature="nl.fitt.gn",
          definition=function(nlfited,control=nlr.control()){
            if(is.Fault(nlfited)) return(nlfited)
            result <- NULL
            v <- attr(nlfited@gpredictor,"gradient")
            hs <- attr(nlfited@gpredictor,"hessian")
            res <-as.numeric(nlfited@gresponse) -as.numeric(nlfited@gpredictor) #residuals(nlfited)
            sigma <- nlfited@scale # old one parameters[["sigma"]]
            grtgr <- t(v) %*% v      						##   F.' * F.     P*P
            ht<-eiginv(grtgr,symmetric =T)						##   (F.' F.)^-1  P*P
            vbar <- apply(v,2,mean)								##   Mean over columns of v, 
            if(class(ht)=="Fault") return(ht)
            covmat <- var(v)										##   (C) Variance covariance matrix.
            cinv <- indifinv(covmat,symmetric =T,stp=F)
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
            if(control$JacobianLeverage == "classic")
              JL <- jaclev(v,hs,res)           ### cllasic for robust use
            else 
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
            ################################################################### output   #########################
            result1 <- structure(.Data=list(
              vmat,
              nl.robmeas(measure=studres,cutofpoint=c(2.5,3,-2.5,-3),name="Studentised Residuals"),
              nl.robmeas(measure=sqrt(elliptnorm),cutofpoint=1,name="Elliptic Norm (Cook Dist)"),
              nl.robmeas(measure=mahd,cutofpoint=c(ctfmahdr,qchisq(.95,p)),name="Regression Mahalanobis Distance"),
              nl.robmeas(measure=mahdata,cutofpoint=ctfmahd1,name="Mahalanobis MVE, data"),
              nl.robmeas(measure=mahdatax,cutofpoint=ctfmahd2,name="Mahalanobis MVE, xs"),
              nl.robmeas(measure=hadi,cutofpoint=c(ctfhadi1,ctfhadi2),name="Hadi potential"),
              nl.robmeas(measure=potmah ,cutofpoint=ctfpotmah,name="Potential mahalanobis"),
              g2,g3
            ),
            .Names=c(
              "vmat",
              "studres",
              "cook",
              "mahd.v",
              "mahd.dt",
              "mahd.xs",
              "hadi",
              "potmah",
              "mvedta","mvex"
            )
            )
            
            if(class(JL)!="Fault"){ 
              result2 <- structure(.Data=list(
                jvmat,
                nl.robmeas(measure=jl.studres,cutofpoint=c(2.5,3,-2.5,-3),name="Jacobian Leverage Studentised Residuals"),
                nl.robmeas(measure=sqrt(jl.elliptnorm),cutofpoint=1,name="Jacobian Leverage Elliptic Norm (Cook Dist)"),
                nl.robmeas(measure=jl.hadi,cutofpoint=c(jl.ctfhadi1,jl.ctfhadi2),name="Jacobian Leverage Hadi potential")
              ),
              .Names=c(
                "jl.vmat",
                "jl.studres",
                "jl.cook",
                "jl.hadi"
              )
              )
              result1[names(result2)] <- result2
            }
            
            return(result1)
          }
)


#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.fitt.gn'                                |
#|                                                                                 |
#|                              Sep 2009, update Nov 2012                          |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
