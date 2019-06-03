# ******************************************** #
# **   Least median square fitt,               #
# ** Optimizing using optimize function        #
# **   that is golden section method.          #
# ******************************************** #
nl.lmsNM<-function(formula,data,start){
  if(formula$p==1){
   
    medsq <- function(par, data){
      datalist<-c(data,start)
      nlmodel <- eval(formula,datalist)
      pred <- as.numeric(nlmodel$predictor)
      resp <- as.numeric(nlmodel$response)
      errorsq <- (resp-pred)^2
      median(errorsq)
    }
    p1 <- unlist(start)
    p2<-p1-.2*p1
    p3<-p1+.2*p1
    fit=optimize( medsq, interval=c(p2,p3),data = data)
    result <- as.list(fit$minimum)
    names(result)<- names(start)
    return(result)
  }
  else{
      medsq <- function(par, data){
        datalist<-c(data,start)
        nlmodel <- eval(formula,datalist)
        pred <- as.numeric(nlmodel$predictor)
        resp <- as.numeric(nlmodel$response)
        errorsq <- (resp-pred)^2
        median(errorsq)
      } 
    fit=optim(par=start,fn= medsq, data =data)
    as.list(fit$par)
  }
}
