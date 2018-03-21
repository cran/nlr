nl.lts <-function(formula,data,start,h=NULL,control=nlr.control()){
  fit1 <- nlsqr(formula,data,start,control)  
  theta1 <- unlist(fit1$parameters[names(formula$par)]) # withough sigma
  eofl <- F
  ri1 <- residuals(fit1)
  n <- length(ri1)
  if(is.null(h)){
    alpha=0.25
    h<-ceiling(n*(1-alpha))  # bigger integer
  }
  dun <- round(runif(h)*n)
  datalist<-NULL
  risq <- ri1^2
  riord <- sort(risq)
  rilow <- risq<=riord[h]
  for(i in 1:length(formula$dependent)) datalist[[formula$dependent[i]]]<-data[[formula$dependent[[i]]]][rilow]
  for(i in 1:length(formula$independent)) datalist[[formula$independent[i]]]<-data[[formula$independent[i]]][rilow]
  fit2 <- nlsqr(formula,datalist,theta1,control)
  plot(fit2)
  
  while(eofl==F){
    risq <- ri1^2
    riord <- sort(risq)
    rilow <- risq<=riord[h]
    for(i in 1:length(formula$dependent)) datalist[[formula$dependent[i]]]<-data[[formula$dependent[[i]]]][rilow]
    for(i in 1:length(formula$independent)) datalist[[formula$independent[i]]]<-data[[formula$independent[i]]][rilow]
    fit2 <- nlsqr(formula,datalist,theta1,control)
    plot(fit2)
    theta2 <- unlist(fit2$parameters[names(formula$par)]) # withough sigma
    diff <- sqrt(sum((theta1-theta2)^2))
    if(diff<=control$tolerance) eofl=T
    theta1 <- theta2
  }
  return(fit2)
}