
#+##################################################################################+
#**    nl.fitt: object. store fitted values of a method.                            \
#**       Slots:                                                                     \
#**             parameters:  list of parameter estimates.                             \
#**             scale:       is sigma, before was in parameters but now added as a    |
#**                          slot. It is standard deviation not variance.             \
#**             correlation: correlation of parameter estimates, computed from         \
#**                          linear approximation.                                      |
#**             formula:     "nl.form" object of called form.                           |
#**             response:    computed response, left hand formuls.                      |
#**             predictor:   computed predictor right hand formula, contains            |
#**                          extra attributes,  gradient and hessian of the            /
#**                          predictor formula.                                       /
#**             curvature:   liniear prediction identifier, list of pe(parameter     |
#**                          effect), IN (Intrinsic nonlinearity), A (A matrix).     |
#**                          see Bates and Watts 1988.                               |
#**                          for robust cases and others it can be meaningful.        \
#**                          its better to define it as an ellement ($) not slot (@)   \
#**             history:     history matrix of itterations, depends on the method.      \
#**                          with named columns.  in modified gauss newton contains,    |
#**                          'iteration', iteration number.                             |
#**                          'ssq', sum of square at each itteration.                   \
#**                          '.....' parameters estiamtes at each column.               \
#**                          'converge' converge factor, defined by batted and Watts.    \
#**                          in definition its NULL in the begining, and diffrent         |
#**                          methods return diffrent (named) column of the matrix,        |
#**                          it can be completlly general.                               /
#**             method:      'fittmethod' object of method used.                        /
#**             Fault:       'list' of Faults in program, depends on methods contains: /
#**                               FL: 'logical', if True (T) error, other Slots will   \
#**                                   be NULL, if False (F) no error method converged   \
#**                               FN: 'integer' error number, if equal(0) zero no error, \ 
#**                                   otherwise is the number of error returned back by   \
#**                                   the method, some times the number shows an warning  |
#**                                   it means output exist i.e (FL=F), but some          |
#**                                   convergence problem occured. the number must be     |
#**                                   defined inside the method.                         /
#**                                                                                     |
########################################################################################+

#+--------------------------------------------------+
#|     class:  nl.fitt                              |
#+--------------------------------------------------+
setClass(Class="nl.fitt",representation=representation(
										parameters =     "list",
										scale =          "numericorNULL",
										correlation =    "numericorNULL",
										form =           "nl.form",
										response =       "vectororMatrix",
										predictor =      "vectororMatrix",
										curvature =      "listorNULL",
										history =        "matrixororNULL",
										method =         "fittmethodorNULL",
										data =           "list",
										sourcefnc =      "callorNULL",
										Fault =          "Fault",
										others =         "listorNULL"
									),
         prototype = prototype(Fault=Fault())
			)
setClassUnion("nl.fittorNULL", c("nl.fitt", "NULL"))

###################################################
##   constructor function, 
##        if length of parameters is zero 
##        it means object is empty or null
###################################################
nl.fitt<-function(parameters ,scale=NULL,correlation=NULL,form ,response =NULL,predictor ,
						curvature =NULL,history =NULL,method=NULL,data,sourcefnc=NULL,Fault,others=NULL){
	if(Fault@FL) return(new("nl.fitt",Fault=Fault))
	if(is.data.frame(history)) history <- as.matrix(history)

  result <- new("nl.fitt",
							parameters =   parameters ,
							scale =        scale,
							correlation=   correlation,
							form =         form, 
							response =     response, 
							predictor =    predictor, 
							curvature =    curvature,
							history =      history,
							method =       method,
							data =         data,
							sourcefnc =    sourcefnc,
							Fault =        Fault,
							others =       others
				)
	

  return(result)
}

###################################################
##   Method $, nl.fitt.
###################################################
#method.skeleton("$","nl.fitt")

setMethod("$","nl.fitt",
	function(x,name){
		slot(x,name)
	} 
)

###################################################
##   Method "$<-", nl.fitt. Its assingment.
###################################################
#method.skeleton("$<-","nl.fitt")  
#method.skeleton("predict","nl.fitt")  
#method.skeleton("residuals","nl.fitt")  
#method.skeleton("recalc","nl.fitt") 
#setMethod("$<-","nl.fitt",
#	function(nl.fitt,name,newvalue){
#		slot(nl.fitt,name) <- newvalue
#	}
#)
###################################################
##   Method residuals, nl.fitt.
###################################################
setMethod(f="residuals",signature="nl.fitt",
	definition=function(object,...){
		dotlist <- list(...)
		if(is.null(dotlist$data)){
			rs <- as.numeric(object@response)
			pr <- as.numeric(object@predictor)
			result <- (rs-pr)
		}
		else{
			datalist <- as.list(dotlist$data)
			datalist[names(object@form$par)] <- object@parameters[names(object@form$par)]
			fmod <- eval(object@form,datalist)
			rs <- as.numeric(fmod@response)
			pr <- as.numeric(fmod@predictor)
			result <- (rs-pr)
		}
	return(result)
	}
)

###################################################
##   Method predict, nl.fitt.
###################################################
setMethod(f="predict",signature="nl.fitt",
	definition=function(object,...){
		dotlist <- list(...)
		if(is.Fault(object)) return(NULL)
		if(is.null(dotlist$newdata)) {
				prnt <- as.list(object$data)
				prnt[names(object@parameters)] <- object@parameters
			  prd <- eval(object@form,prnt)
			  return(prd$predictor)
		}
		else{
			prnt <- as.list(dotlist$newdata)
			prnt[names(object@parameters)] <- object@parameters
      if(is.null(prnt[[object@form@dependent]])) prnt[[object@form@dependent]] <- prnt[[object@form@independent]] # just formality
 			prd <- eval(object@form,prnt)
			return(prd$predictor)
		}
	}
)

###################################################
##   Method hat, nl.fitt.
###################################################
setMethod(f="hat",signature="nl.fitt",#signature(x="nl.fitt",intercept=NULL),
	definition=function(x="nl.fitt",intercept=NULL){
		g <- attr(x$predictor,"gradient")
		gg <- t(g) %*% g
		ggi <- indifinv(gg,stp=F)
		if(is.Fault(ggi)) return(ggi)
		catcher <- ggi %*% t(g)
		hat <- g %*% catcher
		return(hat)
	}
)

###################################################
##   Method plot, nl.fitt.
###################################################
#method.skeleton("plot","nl.fitt") 
setMethod("plot",signature(x="nl.fitt"),
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
    
    if(is.Fault(x)) {
			cat("nl.fitt object returned an error Fault=T   ",str(x$Fault),"\n")
			return()
		}
    if(length(x@form$independent) == 1){
					if(! singlePlot) par(mfrow=c(1,2),...)
					if(is.null(length.out)){
							prdt <- predict(x,newdata=x@data)
							prdt<-as.numeric(prdt)
							prdtx <- x@data[[x@form@independent]]
							orsrt <- order(x@data[[x@form$independent]])
							predint<- predictionI(x)
					} 
					else {
							prdtx <- seq(min(x@data[[x@form@independent]]),max(x@data[[x@form@independent]]),length.out=length.out)
							newdata <- NULL
							newdata[[x@form@independent]] <- prdtx
							prdt <- predict(x,newdata=newdata)
							prdt<-as.numeric(prdt)
							newdata[[x@form@dependent]] <- as.numeric(prdt)							
							orsrt <- order(newdata[[x@form@independent]])
							predint<- predictionI(x,data=newdata)
						}
	
					ymax <- max(x@data[[x@form$dependent]],predint$upperbound,predint$lowerbound)
					ymin <-  min(x@data[[x@form$dependent]],predint$upperbound,predint$lowerbound)
					par(mfrow=c(1,2))
					plot(x@data[[x@form@independent]] , x@data[[x@form@dependent]],
						main=paste(x@form@name,"\n",sub=x@method@method),
						xlab=x@form@independent,
						ylab=x@form@dependent,
						ylim=c(ymin,ymax),...)

					lines(prdtx[orsrt],prdt[orsrt],lty=1)

					lines(prdtx[orsrt],predint$upperbound[orsrt],col=2)
					lines(prdtx[orsrt],predint$lowerbound[orsrt],col=2)
					
					#lines(prdtx[orsrt],predint$upperbound[orsrt],lty=2)
					#lines(prdtx[orsrt],predint$lowerbound[orsrt],lty=2)
					
					plot(x@data[[x@form@independent]] ,residuals(x),main=x@form@name,ylab="Residuals",xlab="Order",...)
					par(mfrow=c(1,1))
					if(history & nrow(x@history) != 0){
						xt <- paste(x@form@independent," nof iteration=",as.character(nrow(x@history))) 
						dev.new()
						par(mfrow=c(2,2))
						plot(x@data[[x@form@independent]] , x@data[[x@form@dependent]],
							main=paste(x@form@name, "\n convergence ,",x@method@method),
							xlab = xt,
							ylab=x@form@dependent,...)
						for(i in 1:nrow(x@history)){
							datalist <- as.list(x@data)
							datalist[ unlist(dimnames(x@history)[2])  ]<- x@history[i,]
							fm <- eval(x@form,datalist)
							lines(x@data[[x@form@independent]][orsrt],as.numeric(fm$predictor)[orsrt])
						}
						plot(x@history[,"objfnc"],type="b",
							main="ssq, sum of square",
							xlab = xt,
							ylab="Sum of Square",...)
	
						plot(x@history[,"converge"],type="b",
							main="converge",
							xlab = xt,
							ylab="Convergance",...)
						par(mfrow=c(1,1))
				
			}
		}
		else{
		      dev.new()
					pairs(data.frame(  x@data[c(x@form@dependent,x@form@independent)]  ))
					plot(residuals(x),main=x@form@name,ylab="Residuals",xlab="Order",...)
					if(history & nrow(x@history) != 0){
							xt <- paste(x@form@independent," nof itteration=",as.character(nrow(x@history))) 
							if(! singlePlot) par(mfrow=c(2,1))
							plot(x@history[,"objfnc"],type="b",
							main="ssq, sum of square",
							xlab = xt,
							ylab=x@form@dependent,...)
							plot(x@history[,"converge"],type="b",
								main="converge",
								xlab = xt,
								ylab=x@form@dependent,...)
							par(mfrow=c(1,1))
			}			
		}
	}
)
#########################################################################
##
##   Method mplot, nl.fitt, plot multiple models
##     Added to this object att 27/11/2012
##   input:
##          mlist: list of object  models
##          case=1, common x
##          case=2, different x
##
########################################################################
"mplot"<-
	function(mlist,case=1,length.out=NULL,...){
			nobj<-length(mlist)
			dot.args <- list(...)
			if(is.null(dot.args$main)) main <- paste(mlist[[1]]@form@name,"\n",sub=mlist[[1]]@method@method)
			else main <- dot.args$main
			if(is.null(dot.args$xlab)) xlab <-mlist[[1]]@form@independent
			else xlab <- dot.args$xlab
			if(is.null(dot.args$ylab)) ylab <- mlist[[1]]@form@dependent
			else ylab <- dot.args$ylab
			
			if(length(mlist[[1]]@form@independent) != 1) return()
			switch(case,
				"1"={
							orsrt <- order(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							xlow<-min(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							ylow<-min(mlist[[1]]@data[[mlist[[1]]@form@dependent]])
							xhigh<-max(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							yhigh<-max(mlist[[1]]@data[[mlist[[1]]@form@dependent]])
							prdt <- as.matrix(as.numeric(predict(mlist[[1]],newdata=mlist[[1]]@data)))
							for(i in 2:nobj){
							  prdt <- cbind(prdt,as.numeric(predict(mlist[[i]],newdata=mlist[[i]]@data)))
							}
							ymax <- max(mlist[[1]]@data[[mlist[[1]]@form@dependent]],prdt)
							ymin <- min(mlist[[1]]@data[[mlist[[1]]@form@dependent]],prdt)
							
							plot(mlist[[1]]@data[[mlist[[1]]@form@independent]] , mlist[[1]]@data[[mlist[[1]]@form@dependent]],
							     main=main,
							     xlab=xlab,
							     ylab=ylab,
							     col=1,
							     ylim=c(ymin,ymax)
							     )
              lines(mlist[[1]]@data[[mlist[[1]]@form@independent]][orsrt],prdt[orsrt,1],col=1,lty=1,...)
							legtxt<-NULL
							legtxt[1]<-mlist[[1]]@method@method
							for(i in 2:nobj){
								lines(mlist[[i]]@data[[mlist[[i]]@form@independent]][orsrt],prdt[orsrt,i],col=1,lty=i,...)
								legtxt[i]<-mlist[[i]]@method@method
							}
              
							legend(xlow,yhigh,legtxt,lty=1:nobj,col=1)#:nobj)
					},
				"2"={
							xlow<-min(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							ylow<-min(mlist[[1]]@data[[mlist[[1]]@form@dependent]])
							xhigh<-max(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							yhigh<-max(mlist[[1]]@data[[mlist[[1]]@form@dependent]])
							for(i in 2:nobj){
								xlow<-min(xlow,mlist[[i]]@data[[mlist[[i]]@form@independent]])
								ylow<-min(ylow,mlist[[i]]@data[[mlist[[i]]@form@dependent]])
								xhigh<-max(xhigh,mlist[[i]]@data[[mlist[[i]]@form@independent]])
								yhigh<-max(yhigh,mlist[[i]]@data[[mlist[[i]]@form@dependent]])
							}
							orsrt <- order(mlist[[1]]@data[[mlist[[1]]@form@independent]])
							plot(mlist[[1]]@data[[mlist[[1]]@form@independent]] , mlist[[1]]@data[[mlist[[1]]@form@dependent]],
								main=main,
								xlab=xlab,
								ylab=ylab,
								xlim=c(xlow,xhigh),ylim=c(ylow,yhigh),col=1,...)
							prdt <- predict(mlist[[1]],newdata=mlist[[1]]@data)
							lines(mlist[[1]]@data[[mlist[[1]]@form@independent]][orsrt],
								prdt[orsrt],col=1,lty=1,pch=1
								)
							legtxt<-NULL
							legtxt[1]<-mlist[[1]]@method@method
							for(i in 2:nobj){
								points(mlist[[i]]@data[[mlist[[i]]@form@independent]] , mlist[[i]]@data[[mlist[[i]]@form@dependent]],col=1,pch=i,...)
								prdt <- predict(mlist[[i]],newdata=mlist[[i]]@data)
								lines(mlist[[1]]@data[[mlist[[i]]@form@independent]][orsrt],prdt[orsrt],col=1,lty=i,...)
								legtxt[i]<-mlist[[i]]@method@method
							}
							legend(xlow,yhigh,legtxt,lty=1:nobj,col=1:nobj,pch=1)#:nobj)
					}
			)
	}

###################################################
#######################################################
##   Method atypicals, nl.fitt, parameter Infrences
##     Added to this object att 25/11/2008
##  According to riaz2009 et al only studres and elliptnorm
##     can detect outliers. 
##     Other measures shown here for feature development.
#######################################################
setGeneric("atypicals", 
           def=function(nlfited,control=nlr.control()) standardGeneric("atypicals"))

setMethod(f="atypicals",signature="nl.fitt",
          definition=function(nlfited,control=nlr.control()){
            if(is.Fault(nlfited)) return(nlfited)
            result <- NULL
            v <- attr(nlfited@predictor,"gradient")
            hs <- attr(nlfited@predictor,"hessian")
            res <- residuals(nlfited)
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

#######################################################
##   Method atypicals.deleted, nl.fitt, parameter Infrences
##     Splited from original atypical function
##      at 17/02/2015, UPM, visitor
##  such methods are not recomended.
##  According to riaz2009 et al only studres and elliptnorm
##     can detect outliers. 
##     These deleted measures shown here for feature development.
#######################################################
setGeneric("atypicals.deleted", 
           def=function(nlfited,control=nlr.control()) standardGeneric("atypicals.deleted"))

setMethod(f="atypicals.deleted",signature="nl.fitt",
          definition=function(nlfited,control=nlr.control()){
            if(is.Fault(nlfited)) return(nlfited)
            result <- NULL
            v <- attr(nlfited@predictor,"gradient")
            res <- residuals(nlfited)
            p <- nlfited$form$p
            n <- length(res)
            sigma <- nlfited@scale # old one parameters[["sigma"]]
            grtgr <- t(v) %*% v      						##   F.' * F.     P*P
            ht<-eiginv(grtgr,symmetric =T)						##   (F.' F.)^-1  P*P
            dataset <- nlfited$data[c(nlfited$form$independent,nlfited$form$dependent)]	
            dataset <- data.frame(dataset)
            datasetx <- nlfited$data[c(nlfited$form$independent)]
            datasetx <- data.frame(datasetx)
            ##  compute vmat: square of mahalanobis distance
            vmat <- rep(0,n)
            for(i in 1:nrow(v)){									##   Cllasic, center and variance are classic.
              b1<-v[i,]											##  vi', will store in row, 1*p
              b2<-b1 %*% ht										##  vi' (g'g)^-1            1*p
              vmat[i] <- b2 %*% b1								##  wii = vi' (g'g)^-1 vi   1*1
            }
            
            #####################################
            if(control$JacobianLeverage == "classic")
              JL <- jaclev(v,attr(nlfited@predictor,"hessian"),res)           ### cllasic for robust use
            else 
              JL <- JacobianLeverage(nlfited)					###   Jacobian Leverage
            if(class(JL)!="Fault"){
              jvmat <- diag(JL)
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
              newfiti <- recalc(nlfited,data=dst2)					##   fit all (-i)
              newfiti@Fault@FF <- "nlfited.atypicals"
              if(newfiti$Fault$FL) return(newfiti)
              d.sigma[i]  <- newfiti$scale 	##  sigma(-i)
              yhati <- predict(newfiti,newdata=dataset[i,])				##  yhat(-i) grd&hess
              d.yhat[i] <- as.numeric(yhati)						##  yhat(-i)
              dfbt <- unlist(nlfited$parameters[names(nlfited$parameters)!="sigma"]) - 
                unlist(newfiti$parameters[names(newfiti$parameters)!="sigma"])
              d.fbetas[i,] <- dfbt / newfiti$scale / sqrt(diag(ht))   ##   DFBETASi = bj-bj(i) / sg(i) sart(x'x-1 jj)
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
              d.yhat,
              nl.robmeas(measure=delstud ,cutofpoint=c(2.5,3,-2.5,-3),name="Deletion Studentized"),
              nl.robmeas(measure=d.ffits ,cutofpoint=ctfdffits ,name="DFFITS"),
              nl.robmeas(measure=atkinson ,cutofpoint=c(2,ctfatk),name="Atkinson Distance"),
              d.fbetas
            ),
            .Names=c(
              "yihat",
              "delstud",
              "dffits",
              "atk",
              "dfbetas")
            )
            
            if(class(JL)!="Fault"){ 
              result2 <- structure(.Data=list(
                nl.robmeas(measure=jl.delstud ,cutofpoint=c(2.5,3,-2.5,-3),name="Jacobian Leverage Deletion Studentized"),
                nl.robmeas(measure=jl.d.ffits ,cutofpoint=ctfdffits ,name="Jacobian Leverage DFFITSi"),
                nl.robmeas(measure=jl.atkinson ,cutofpoint=c(2,jl.ctfatk),name="Jacobian Leverage Atkinson Distance")
              ),
              .Names=c(
                "jl.delstud",
                "jl.dffits",
                "jl.atk")
              )
              result1[names(result2)] <- result2
            }
            
            return(result1)
            
          }
)


#######################################################
##   Method recalc, nl.fitt, recalculate object "object
##   for new values "colin".  added 28/12/2008
#######################################################
setGeneric("recalc", 
	def=function(object,...) standardGeneric("recalc"))
setMethod(f="recalc",signature="nl.fitt",
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
##   Method JacobianLeverage, nl.fitt, parameter Infrences
##     Added to this object att 23/3/2009
#######################################################
setGeneric("JacobianLeverage", 
	def=function(nlfited) standardGeneric("JacobianLeverage"))

setMethod("JacobianLeverage",
          signature(nlfited = "nl.fitt"),
          function (nlfited) 
           {
		          if(is.Fault(nlfited)) return(nlfited)
		          gradient<-attr(nlfited@predictor,"gradient")
              hessian<-attr(nlfited@predictor,"hessian")
		          rsd<-residuals(nlfited)
              
		          prva <-   prodVA(hessian,rsd)                ##  [e] [w]   p*p
		          .expr1 <- t(gradient) %*% gradient           ##  vT * v    p*p
		          .expr2 <- .expr1 - prva                      ##  vT v - [e] [w]    p*p
		          .expr3 <- indifinv(.expr2,stp=F,symmetric=T)   ##  (vT v - [e] [w])^-1   p*p
		          if(class(.expr3)=="Fault") return(.expr3)
		          .temp1 <- gradient %*% .expr3                ##  v * (vT v - [e] [w])^-1    n*p
		          JacLeverage <- .temp1 %*% t(gradient)        ##  v * (vT v - [e] [w])^-1 vT  n*n
		          return(JacLeverage)
		          
	          }
)


###################################################
##   Method parInfer, nl.fitt, parameter Infrences
###################################################

setGeneric(name="parInfer",def=function(object,confidence=.95) standardGeneric("parInfer"))

setMethod("parInfer",signature(object = "nl.fitt"),
  function (object, confidence=0.95) 
  {
			
      .expr1 <- attr(object@predictor,"gradient")							##  F
			.expr2 <- t(.expr1) %*% .expr1										##  F' F

			.expr3 <- eiginv( .expr2 ,stp=F)
			if(is.Fault(.expr3)) .expr3 <- ginv(.expr2)
			covmatrix <- object@scale ^ 2 * .expr3 	##  sgm2 * (F'F)^-1 
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
setGeneric("predictionI", 
	def=function(nlfited,confidence=.95,data=NULL) standardGeneric("predictionI"))

setMethod("predictionI",signature(nlfited="nl.fitt"),
	function(nlfited,confidence=.95,data=NULL){
		if(is.null(data)){
			yhat <- predict(nlfited)
			grdy <- attr(yhat,"gradient")
			yhat <- as.numeric(yhat)
			m <- nrow(grdy)
			data <-nlfited@data[[nlfited@form$independent]]
		}
		else{
			yhat <- predict(nlfited,newdata=data)
			grdy <- attr(yhat,"gradient")
			yhat <- as.numeric(yhat)			
			m <- nrow(grdy)
		}
		sigma <- nlfited$scale
		gpg<-parInfer(nlfited)$covmat             ## sigma embded here
		
		t6<-rep(0,times=m)
		for(i in 1:m){
			t5<-grdy[i,]%*%gpg			#*** 1*p
			t6[i]<-t5%*%(grdy[i,])       #*** 1*1
		}
		rf <- qt((1+confidence)/2,m-nlfited@form$p)    #### note: calculate upper = 1-alpha/2 = (1+conf)/2 .... qt calculate lower bound
		vp <- sigma^2+t6
		svp <- sqrt(vp)
		up<-yhat+rf*svp
		lw<-yhat-rf*svp

		return(list(lowerbound=lw,upperbound=up))
	}
)



#######################################################+
##   Method ACF, nl.fitt, Empirical auto correlation  /
##   function of residuals                            |
##     Added to this object att 09/07/2014             \
##   it is s3 function style, same as ACF.lme           |
#######################################################-+

#method.skeleton("acf", "nl.fitt")
setMethod("acf",
          signature(x = "nl.fitt"),
          function (x, lag.max = NULL, type = c("correlation", "covariance", 
                                                "partial")[1], plot = TRUE, na.action = na.fail, demean = TRUE, 
                    drop.lag.0 = TRUE,...)
          {
            rsd <- residuals(x)
            if(is.null(lag.max)) lag.max <- length(rsd)
            return(acf(rsd),type=type,plot=plot,na.action=na.action,demean=demean,...)
          }
)
#############################################################|
##   Method replacement, nl.fitt.                             \

###################################################            \
#--------------------------------------------------------------/
setReplaceMethod(f="[", signature="nl.fitt",
                 definition=function (x,i,value){
                   if (i=="others"){x@selfStart<-value}
                   validObject(x)
                   return(x)
                 }  
)
#same--------------------------------------------------------------/
setReplaceMethod(f="$", signature="nl.form",
                 definition=function (x,name,value){
                   if (name=="selfStart"){x@selfStart<-value}
                   else stop("There is no slot called",name)
                   return(x)
                 }  
)
#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.fitt'                                   |
#|                                                                                 |
#|                                                                                 |
#|                            09/2008                                              |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
