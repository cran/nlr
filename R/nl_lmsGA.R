nl.lmsGA<-function(formula,data,start,min=NULL,max=NULL,type="real-valued"){
  
  medsq <- function(par){
    names(par) <- names(formula@par)
    datalist<-c(par,data)
    nlmodel <- eval(formula,datalist)
    if (is.Fault(nlmodel)) return(NA)
    pred <- as.numeric(nlmodel$predictor)
    resp <- as.numeric(nlmodel$response)
    errorsq <- (resp-pred)^2
    -median(errorsq)
  } 
  if ( is.null(min)){
    .temp <- unlist(start)
    min<-.temp*.9
    max<-.temp*1.1
  }
  fitGA <- ga(type = type, fitness = function(par) medsq (par), min = min, max = max,popSize = 50,keepBest=TRUE, 
              monitor = FALSE)
  MedianSquares <- as.numeric(NULL)  
  
  i<-1
  par=as.list(summary(fitGA)$solution[i,])
  names(par)<-names(formula$par)
  dataset<-c(par,data)
  y2<-eval(formula,dataset)
  y3<-as.numeric(y2$predictor)
  y4 <- as.numeric(y2$response)
  ri <- y3-y4
  MedianSquares[i]<-(median(ri^2))
  parameters <- par
  rilms<-ri
  response <- y2$response
  predictor<- y2$predictor
  htheta <- MedianSquares[i]
  if(nrow(summary(fitGA)$solution) >1){
    for(i in 2:nrow(summary(fitGA)$solution)){
      par=as.list(summary(fitGA)$solution[i,])
      names(par)<-names(formula$par)
      dataset<-c(par,data)
      y2<-eval(formula,dataset)
      y3<-as.numeric(y2$predictor)
      y4 <- as.numeric(y2$response)
      ri <- y3-y4
      MedianSquares[i]<-(median(ri^2))
      if (MedianSquares[i] <= fitGA@fitnessValue) {
        parameters <- par
        rilms<-ri
        response <- y2$response
        predictor<- y2$predictor
        htheta <- MedianSquares[i]
      }
    }
  }
  n <-length(ri)
  history <- summary(fitGA)$solution
  dimnames(history)[[2]] <- names(formula$par)
  names(MedianSquares) <- "objfnc"
  history <- cbind(history,MedianSquares)
  
  history <- cbind(1:nrow(history),history)
  dimnames(history)[[2]][1]<-"iteration"
  
  # calculate scale 
  s0 <- 1.4826*(1+(5/(length(n-formula@p)))) * sqrt(rilms^2)
  ri2 <- rilms[rilms<=2.5*s0]
  n2 <- length(ri2)
  sigma2 <- sum(ri2^2)/(n2-formula@p)
  sigma <- sqrt(sigma2)
  ybar <- mean(response)
  nlrho <- 1 - sum( ( response - predictor )  ^ 2 ) / 
    sum( (predictor - ybar ) ^ 2 )
  curv1 <- curvature(gradient = attr(predictor,"gradient"),
                     hessian = attr(predictor,"hessian"),
                     sigma = sigma)
  
  result <- nl.fitt.rob(parameters =  parameters,
                        scale =        sigma,
                        correlation =  nlrho,
                        form =         formula,
                        response =     response,
                        predictor =    predictor,
                        curvature =    curv1,
                        history =      history,
                        method =   	 fittmethod(methodID=			12,
                                                methodBR=			16,       ### Genetic algorithm
                                                detailBR=			"Genetic Algorithm",
                                                subroutine=		"nl.lmsGA"),
                        data =         as.list(data),
                        sourcefnc =     match.call(),
                        Fault =        Fault(),
                        htheta =       htheta,     # objective function
                        rho =          rilms^2,      # loss values
                        ri =           rilms,
                        curvrob =      NULL,
                        others =       list(fitGA=fitGA)
  )
  return(result)
}


