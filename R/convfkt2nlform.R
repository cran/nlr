convfkt2nlform<-function(fktlistex,namesdata=NULL){
  fktmodels <- fktlistex
  a0 <- call("~", fktmodels$fkt)
  a1 <- formula(a0)
  a2 <- names(fktmodels$par)
  b <- deriv3(a1, a2)
  a5 <- call("~", (~yr)[[2]], b)
  fktmodels[[9]] <- a5
  names(fktmodels)[9] <- "fkm"
  .org1 <- call("~",(~yr)[[2]],a0[[2]])
  fktmodels[[10]] <- .org1
  names(fktmodels)[10] <- "origin"		
  result <- nl.form(form=fktmodels$fkm,p=length(fktmodels$par),
        inv=fktmodels$inv,par=fktmodels$par,name=fktmodels$name,
        independent="xr",dependent="yr",origin=fktmodels$origin)
  return(result)
}

