#+#################################################################################+
#**       nl.fitt.rob:  a nonlinear fitted robust method. for storing fitted vlues  |
#**				Inherit from nl.fitt with all Methods.                                     \
#**			Extra Slots                                                                  /
#**             htheta:  h(t) = sum rh( ri / sigma), loss rho function.             /
#**             rho:     ho function, rho(ri / sigma),                             /
#**             
#+--------------------------------------------------+
#|      class:  nl.fitt.rob                         |
#+--------------------------------------------------+
setClass("nl.fitt.rob",representation(
									"nl.fitt",
									htheta =      "vectororNULL",
									rho =         "vectororNULL",
									ri =          "vectororNULL",
 									curvrob =     "listorNULL",
									robform =     "nl.formorNULL"
									)
		)
setClassUnion("nl.fitt.roborNULL", c("nl.fitt.rob", "NULL"))
###################################################
##   constructor
###################################################
nl.fitt.rob<-function(parameters=NULL, scale=NULL, correlation=NULL,form=NULL ,response =NULL,predictor=NULL ,
						curvature =NULL,history =NULL,method=NULL,data=NULL,sourcefnc=NULL,Fault=new("Fault",FL=FALSE,FN=0,FT="",FF=""),others=NULL,
						htheta=NULL,rho=NULL,ri=NULL,curvrob = NULL ,robform=NULL){
  if(Fault@FL) return(Fault) #new("nl.fitt.rob",form=NULL,Fault=Fault))
  
  result <- new("nl.fitt.rob",
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
							htheta =       htheta,
							rho =          rho,
							ri =           ri,
							curvrob =      curvrob,
							robform =      robform
				)
	return(result)
}
###################################################
##   Method plot, nl.fitt.rob
###################################################

setMethod("plot",signature(x="nl.fitt.rob"),
          function(x,y = "missing",control=nlr.control(history=FALSE,length.out=NULL,singlePlot=FALSE),...){
            dotarg <- list(...)
            if (is.null(dotarg$control) )
            {
              control <- nlr.control()
            }
            else
            {
              control=dotarg$control  
            }
            history=control$history
            length.out=control$length.out
            singlePlot=control$singlePlot
            if(x$Fault$FL == T) {
              cat("nl.fitt object returned an error Fault=T   ",str(x$Fault),"\n")
              return()
            }
            
            if(length(x@form@independent) == 1){   ######    univariate #######
                                                   if(x@method@methodID == 12){   # ***** lms, least median square, estimate *****
                                                     plotlms(x,history=history,...)
                                                   }                            #  ***********************************************
                                                   else{
                                                   if(is.null(length.out)){
                                                     prdt <- predict(x,newdata=x@data)
                                                     prdt<-as.numeric(prdt)
                                                     prdtx  <- x@data[[x@form@independent]]
                                                     orsrt <- order(x@data[[x@form@independent]])
                                                     predint<- predictionI(x)
                                                   } 
                                                   else {
                                                     prdtx  <- seq(min(x@data[[x@form@independent]]),max(x@data[[x@form@independent]]),length.out=length.out)
                                                     newdata <- NULL
                                                     newdata[[x@form@independent]] <- prdtx 
                                                     prdt <- predict(x,newdata=newdata)
                                                     newdata[[x@form@dependent]] <- as.numeric(prdt)
                                                     orsrt <- order(newdata[[x@form@independent]])
                                                     prdt<-as.numeric(prdt)
                                                     predint<- predictionI(x,data=newdata)
                                                   }
                                                   
                                                   
                                                   ymax <- max(x@data[[x@form$dependent]],predint$upperbound,predint$lowerbound,
                                                               prdt)
                                                   ymin <-  min(x@data[[x@form$dependent]],predint$upperbound,predint$lowerbound,
                                                                prdt)
                                                   
                                                   if(! singlePlot) par(mfrow=c(1,2),...)
                                                   par(mfrow=c(1,2))
                                                   plot(x@data[[x@form@independent]] , x@data[[x@form@dependent]],
                                                        main=paste(x@form@name,"\n",sub=x@method@method),
                                                        xlab=x@form@independent,
                                                        ylab=x@form@dependent,
                                                        ylim=c(ymin,ymax),...)
                                                   
                                                   lines(prdtx[orsrt],prdt[orsrt],lty=1)
                                                   
                                                   lines(prdtx[orsrt],predint$upperbound[orsrt],col=2)
                                                   lines(prdtx[orsrt],predint$lowerbound[orsrt],col=2)
                                                   
                                                   #lines(prdtx[orsrt],predint$upperbound[orsrt], lty=2)
                                                   #lines(prdtx[orsrt],predint$lowerbound[orsrt],lty=2)
                                                   
                                                   plot(x@data[[x@form@independent]],residuals(x),
                                                        main=x@form@name,
                                                        ylab="Residuals",xlab="Order",...)
                                                   par(mfrow=c(1,1))
                                                   if(history & nrow(x@history) != 0){
                                                     dev.new()
                                                     par(mfrow=c(2,2))
                                                     xt <- paste(x@form@independent," nof iteration=",as.character(nrow(x@history))) 
                                                     mt <- paste(x@form@name, "\n convergence ,",x@method@method)
                                                     #*****************************
                                                     #** plto converge history   **
                                                     #*****************************
                                                     plot(x@data[[x@form@independent]] , x@data[[x@form@dependent]],
                                                          main=mt,
                                                          xlab = xt,
                                                          ylab=x@form@dependent,...)
                                                     for(i in 1:nrow(x@history)){
                                                       datalist <- as.list(x@data)
                                                       datalist[ unlist(dimnames(x@history)[[2]])  ]<- x@history[i,]
                                                       fm <- eval(x@form,datalist)
                                                       lines(x@data[[x@form@independent]][orsrt],as.numeric(fm$predictor)[orsrt])
                                                     }
                                                     #*****************************
                                                     #** plot rho history        **
                                                     #*****************************
                                                     ri <- as.numeric(residuals(x))
                                                     ord <- order(ri)
                                                     for(i in 1:nrow(x@history)){
                                                       datalist<-list(NULL)
                                                       datalist<- as.list(x@history[i,])#names(x@form@par)])
                                                       fm <- robloss(x@form,data= x@data,start=datalist,robfunc=x$robform)
                                                       if(i == 1) plot(ri[ord],as.numeric(fm$rho)[ord],
                                                                       main=mt,
                                                                       xlab = xt,
                                                                       ylab="Rho function",...)
                                                       else
                                                         lines(ri[ord],as.numeric(fm$rho)[ord])
                                                     }
                                                     #*****************************
                                                     #** plot objective function **
                                                     #*****************************                                                   
                                                     plot(x@history[,"iteration"],x@history[,"objfnc"],
                                                          xlab="Iteration",ylab="robust loss",type="b")
                                                     #*****************************
                                                     #** plot convergence history**
                                                     #*****************************                                                     
                                                     plot(x@history[,"iteration"],x@history[,"converge"],
                                                          xlab="Iteration",ylab="Convergence",type="b")
                                                     par(mfrow=c(1,1))
                                                   }
                                                   }   ##  end os else, other of lms ##
            }  
            else{    ########  multivariate   #######
                     pairs(data.frame(  x@data[c(x@form@dependent,x@form@independent)]  ))
                     plot(residuals(x),main=x@form@name,ylab="Residuals",xlab="Order",...)
                     if(history & nrow(x@history)){
                       if(! singlePlot) par(mfrow=c(2,2))
                       #*****************************
                       #** plot rho history        **
                       #*****************************
                       xt <- paste(" nof itteration=",as.character(nrow(x@history))) 
                       mt <- paste(x@form@name, "\n convergence ,",x@method@method)  			
                       for(i in 1:nrow(x@history)){
                         dataist<-list(NULL)
                         datalist<- as.list(x@history[i,names(x@parameters)])
                         fm <- robloss(x@form,data= x@data,start=datalist,robfunc=x$robform)
                         ri <- as.numeric(x$ri)
                         ord <- order(ri)
                         if(i == 1) plot(ri[ord],as.numeric(fm$rho)[ord],
                                         main=paste(x@form@name, "\n convergence ,",x@method@method),
                                         xlab = xt,
                                         ylab=x@form@dependent,...)
                         else
                           lines(ri[ord],as.numeric(fm$rho)[ord])
                       }
                       plot(x@history[,"iteration"],x@history[,"objfnc"],
                            xlab="Iteration",ylab="robust loss",type="b")
                       plot(x@history[,"iteration"],x@history[,"converge"],
                            xlab="Iteration",ylab="Convergence",type="b")
                       par(mfrow=c(1,1))
                     }
            }
          }
)


# ######################################################
# #   Method JacobianLeverage, nl.fitt.rob, parameter Infrences
# #     Added to this object att 09/05/2009
# ######################################################
setMethod("JacobianLeverage","nl.fitt.rob",
	function(nlfited){
		if(class(nlfited)=="Fault") return(nlfited)
		if(nlfited$Fault$FL==T) return(nlfited$Fault)
		zeta <- attr(nlfited@rho,"hessian")
		phsi <- attr(nlfited@rho,"gradient")
		vmat <- attr(nlfited@predictor,"gradient")

		warray <- attr(nlfited@predictor,"hessian")
		sigma <- nlfited@scale

		.zhv <- zeta * vmat                   ##  z O V           (n.p)
		.zhvv <- t(.zhv) %*% vmat             ##  (z O V)T V      (p.p)
		.expr1 <- .zhvv / sigma               ##  (z O v)T V / sg (p.p)
		.sw <- prodVA(warray,phsi)            ##  [si] [w]        (p.p)
		.Amat <- .expr1 - .sw                 ##  A matrix = 1/sg (zv)T v /sg -[s][w]
		.Ainv <- indifinv(.Amat,stp=F)        ##  A inverse
		if(class(.Ainv)=="Fault") return(.Ainv)
		.expr2 <- .Ainv %*% t(.zhv)
		JL <- vmat %*% .expr2
		JL <- JL / sigma
		#JL <- jaclev(attr(nlfited@predictor,"gradient"),attr(nlfited@predictor,"hessian"),residuals(nlfited))
		return(JL)
		
	}
)

###################################################
##   Method parInfer, nl.fitt.rob, parameter Infrences
##      upgrade 5/12/2012 Stockholm
###################################################


setMethod("parInfer","nl.fitt.rob",
	function(object,confidence=.95){
		if(is.null(object@scale))
		{
      Fault2 <- Fault(FL=T,FN=8,FT=,
			nlr::db.Fault[nlr::db.Fault$FN==8,]$FT,FF="parInfer")
			result <- Fault2
			return(result)
		}
			.expr1 <- attr(object@predictor,"gradient")							##  F
			.expr2 <- t(.expr1) %*% .expr1											##  F' F
			.expr3 <- sum(attr(object@rho,"gradient"))
			n <- nrow(.expr1)
			p <- object$form$p
			spsi <- sum(.expr3^2) 
			.expr4 <- eiginv( .expr2 ,symmetric=T,stp=F)
			if(is.Fault(.expr4)) .expr4 <- ginv(.expr2)
			spsid <- sum(attr(object@rho,"hessian"))
			covmatrix <- object$scale^2 * ( spsi / spsid^2 ) * (n*n/(n-p)) * .expr4
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
##     Added to this object att 05/12/2012
#######################################################


setMethod("predictionI","nl.fitt.rob",
	function(nlfited,confidence=.95,data=NULL){
		if(is.null(data)){
			yhat <- predict(nlfited)
			data <-nlfited@data[[nlfited@form$independent]]
			grdy <- attr(yhat,"gradient")
			yhat <- as.numeric(yhat)			
			m <- nrow(grdy)
		}
		else{
			yhat <- predict(nlfited,newdata=data)
			grdy <- attr(yhat,"gradient")
			m <- nrow(grdy)
		}
		sigma <- nlfited$scale
		gpg<-parInfer(nlfited)$covmat             ## sigma embded here
    gpgi <- indifinv(gpg,stp=F)

		if(is.Fault(gpgi)) gpgi <- ginv(gpg)
		t6<-rep(0,times=m)
		for(i in 1:m){
			t5<-grdy[i,]%*%gpg			#*** 1*p
			t6[i]<-t5%*%(grdy[i,])       #*** 1*1
		}
		rf <- qt((1+confidence)/2,m-nlfited@form$p)
		vp <- sigma^2+t6
		svp <- sqrt(vp)
		up<-yhat+rf*svp
		lw<-yhat-rf*svp
		return(list(lowerbound=lw,upperbound=up))
	}
)

#######################################################
##   Method recalc, nl.fitt, recalculate object "object
##   for new values "object".  added 17/02/2015
#######################################################
setMethod(f="recalc",signature="nl.fitt.rob",
          definition=function(object="nl.fitt",...){
            dotlist <- list(...)
            objcall <- object@sourcefnc
            if(is.null(objcall)) return(Fault(FN=13,FF="recalc nl.fitt"))
            objcall <- object@sourcefnc
            if(is.null(objcall)) return(Fault(FN=13,FF="recalc nl.fitt"))
            if(is.null(dotlist$data)){
              objcall$data <- quote(object@data)
            }
            else{
              objcall$data<-dotlist$data
            }
            objcall$start <- quote(object$parameters)
            
            renew <- eval(objcall )      ######### wrong need data corect it later########
            return(renew)	
          }
)

#######################################################
##   Method dlev, nl.fitt
##     Added to this object att 17/6/2015
##   calculate dlev (difference in leverage)
#######################################################
setGeneric("dlev", 
           def=function(nlfited) standardGeneric("dlev"))

setMethod("dlev",
          signature(nlfited = "nl.fitt.rob"),
          function (nlfited) 
          {
            if(is.Fault(nlfited)) return(nlfited)
            hatmat <- hat(nlfited)
            jaclevmat <- JacobianLeverage(nlfited)
            
            ad1 <- diag(hatmat)
            ad2 <- diag(jaclevmat)
            dlevdiff <- ad2-ad1
            dlevrelat<- dlevdiff 
            cutofpoint<-median(dlevdiff) - 3 * mad(dlevdiff)
            result <- nl.robmeas(measure=dlevrelat ,cutofpoint=cutofpoint,name="Difference in Leverage DLEV")
            
            return(result)
          }
)

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.fitt.gn'                                |
#|                                                                                 |
#|                              Sep 2009 , modified Nov 2012                       |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
