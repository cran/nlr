 var1<-
function(x)
{
	nn <- length(x)
	sum((x - sum(x)/nn)^2)/(nn - 1)
}
