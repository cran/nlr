#************************************************************************
#**   robust loss unction                                              **
#**     start: argument must be a list containing "sigma".             **
#**     argument ... is the arguments for 'robfunc' extra argument     **
#**     for 'robfunc' the robust rho function.                         **
#**     Fault: = 9,  sigma is not a part of start list.                **
#************************************************************************
robloss.gn <- function(formula, data, start, robfunc, rmat, control=nlr.control(robscale=T),
		...)
{
	Fault2 <- Fault()
	data <- as.list(data)
	start <- as.list(start)
	.tmp <- (is.null(start$sigma))
	if (.tmp)
	{
		return(list(Fault=Fault(FL=T, FN=8, FF="robloss.gn")))
	}
	datalist <- c(start, data)
	fmod <- fmod2 <- eval(formula, datalist)						#  f(t) n*1									#
	if(is.Fault(fmod)){
		result <- fmod2
		result@FF<-"robloss.gn"
		return(result)
	}
	####################### GN
	fmod$predictor <- transformNR(fmod$predictor,rmat)
	fmod$response <- transformNR(fmod$response,rmat)
	###########################	
	prd <- as.numeric(fmod$predictor)						#  f(t) values								#
	n <- length(prd)
	p <- formula$p
	gprd <- attr(fmod$predictor, "gradient")			#  f'(t)	Gradient	n*p						#
	hprd <- attr(fmod$predictor, "hessian")				#  f"(t)	Hessian	n*p*p					#
	rsd <- as.numeric(fmod$response) - prd				#  ri		Residals	n*1						#
	ris <- rsd / start$sigma
	if(control$robscale)	robvalue <- robfunc$fnc(ris/robfunc$arguments$k1, ...)  #  rho(ri/sigma)    n*1  #
	else robvalue <- robfunc$fnc(ris, ...)
	rho <- as.numeric(robvalue)
	htheta <- sum(rho)
	grho <- attr(robvalue, "gradient")						#  V=rho'(ri/tht)		n*1						#
	hrho <- attr(robvalue, "hessian")						#  rho"(ri/tht)			n*1					#
	gres <- -gprd												#  J=r'i=-grad(f)		n*p						#
	hres <- -hprd												#  r"i					n*p*p					#
	temp <- t(gres)
	gradh <- temp %*% grho									#  J'V, gradient of rho loss p*p			#
	gradh <- gradh / start$sigma							#  g(h) = 1/sgm * J'* V						#
	dt <- diag(hrho)											#  D(thta)=diag[  rho"(ri/sgm)  ]   n*1	#
	hessh.p1 <- t(gres) %*% dt								#  J'D											#
	hessh.p1 <- hessh.p1 %*% gres							#  J'DJ										#
	hessh.p1 <- hessh.p1 / (start$sigma ^ 2)				#  hessian part 1 = 1/sigma^2  *  J'*D*J	#
	atheta <- matrix(rep(0, p * p), nrow=p)
	for (i in 1:n)
	{
		a <- grho[i] * hres[i, , ]
		atheta <- a + atheta
	}
	
	atheta <- atheta / start$sigma
	hessh.p2 <- atheta
	hessh <- hessh.p1 + hessh.p2
#	dtilda <- diag(start$sigma * grho / rsd)				#  Dtilda(hteta) = rho' / (ri/sigma)		#
	attr(htheta, "gradient") <- gradh
	attr(htheta, "hessian") <- hessh
	attr(rsd, "gradient") <- gres
	attr(rsd, "hessian") <- hres
	result <- list(htheta=htheta, rho=robvalue, ri=rsd,
			 hessh.p1=hessh.p1, hessh.p2=hessh.p2,
			 fmod=fmod2, Fault=Fault2)
#	result <- list(htheta=htheta, rho=robvalue, ri=rsd,
#			 hessh.p1=hessh.p1, hessh.p2=hessh.p2, dtilda=dtilda,
#			 fmod=fmod2, Fault=Fault2)
	
	#****************************************************
	#*****   output:                                ***** 
	#*****  htheta=sum(rho)  fixed number           *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  rho(ri/sigma)       n*1 vector          *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  ri               residuals              *****
	#*****    contains "gradient" & "hessian" attr  *****
	#*****  h(theta)=he	ss.... p1+p2 hessian of (h)*****
	#*****  hess.p1 = p1                            *****
	#*****  hess.p2=  p2                            *****
	#*****  dtilda, new D(thilda) function          *****
	#*****  fmod:  computed function contains       *****
	#*****  response and or its gradient and hessian*****
	#*****  predictor and or its gradient & hessian *****
	#****************************************************
	return(result)
}
