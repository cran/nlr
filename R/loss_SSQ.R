loss.SSQ<-function(formula, data, start,vm=NULL,rm=NULL,...)
{
	Fault2 <- Fault()
	datalist <- c(data, start)
	datalist = as.list(datalist)
	fmod <- eval(formula, datalist)
	if(is.Fault(fmod)){
	  fmod@FF = "loss.SSQ"
	  return(fmod)
	}
	if(! is.null(vm)) {
		rm <- eiginv(t(chol(vm)))
		fmod$predictor<- transformNR(fmod$predictor,rm)
		fmod$response<- transformNR(fmod$response,rm)
	}
	resp <- as.numeric(fmod$response)
	pred <- as.numeric(fmod$predictor)
  ri <- resp - pred
	value <- sum(ri^2)
  value <- as.numeric(value)
	ybar <- mean(resp)
	correlation <- 1 - value / sum(  ( pred - ybar )^2   )
	result <- list(value = value,correlation=correlation,fmod=fmod)
	return(result)
}
