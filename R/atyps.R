atyps <- function(nlfited){
		v <- attr(nlfited@predictor,"gradient")
		hs <- attr(nlfited@predictor,"hessian")
		res <- residuals(nlfited)
		outliers1 <- nlout(nlfited)
		jlevmat <- jaclev(v,hs,res)
		jlev <- diag(jlevmat)
		r1 <- list(jlev=jlev)
		r2 <- c(r1,outliers1)
		result <- as.list(r2)
		return(result)
}