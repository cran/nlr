
zvalues<-
function(res, ni, xo)
{
	# Differences of 2
	n <- length(res)
	dres <- diff(res)
	z <- c((dres * dres), (dres[n - 1] * dres[n - 1]))/2		
	for(i in (1:length(ni))[ni > 1]){
		z[xo == i] <- var(res[xo == i])
	}
	z
}
#################
rzvalues <- function(res, ni, xo)
{
	# Differences of 2
	n <- length(res)
	dres <- diff(res)
	z <- c((dres * dres), (dres[n - 1] * dres[n - 1]))/2
	for(i in (1:length(ni))[ni > 1]){
		z[xo == i] <- (mad(res[xo == i]))^2
	}
	z
}
