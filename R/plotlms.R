# ****** plot for lms, will be called from nl.fitt.rob.plot function ********
plotlms<-function(x,history=F,length.out=NULL,...){
  if(is.null(length.out)){
    prdt <- predict(x,newdata=x@data)
    prdt<-as.numeric(prdt)
    prdtx  <- x@data[[x@form@independent]]
    orsrt <- order(x@data[[x@form@independent]])
  } 
  else {
    prdtx  <- seq(min(x@data[[x@form@independent]]),max(x@data[[x@form@independent]]),length.out=length.out)
    newdata <- NULL
    newdata[[x@form@independent]] <- prdtx 
    prdt <- predict(x,newdata=newdata)
    newdata[[x@form@dependent]] <- as.numeric(prdt)
    orsrt <- order(newdata[[x@form@independent]])
    prdt<-as.numeric(prdt)
  }
  ymax <- max(x@data[[x@form$dependent]],
              prdt)
  ymin <-  min(x@data[[x@form$dependent]],
               prdt)
  
  par(mfrow=c(1,2),...)
  plot(x@data[[x@form@independent]] , x@data[[x@form@dependent]],
       main=paste(x@form@name,"\n",sub=x@method@method),
       xlab=x@form@independent,
       ylab=x@form@dependent,
       ylim=c(ymin,ymax),
       ...)
  
  lines(prdtx[orsrt],prdt[orsrt],lty=1)
  plot(x@data[[x@form@independent]],residuals(x),
       main=x@form@name,
       ylab="residuals",xlab="order",...)
  par(mfrow=c(1,1))
  if(history & nrow(x@history) != 0){
    par(mfrow=c(1,2))
    xt <- paste(x@form@independent," nof itteration=",as.character(nrow(x@history))) 
    mt <- paste(x@form@name, "\n GA POpulation ,",x@method@method)
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
    #** plot objective function **
    #*****************************                                                   
    plot(x@history[,"iteration"],x@history[,"objfnc"],main="Population MedianSquares",
         xlab="Iteration",ylab="robust loss",type="b")
    
    par(mfrow=c(1,1))
  }
}