#* +-------------------------------------------------------------------------------+ *
#* |    parInfer: Extend cllasical linear regression parameter inference           | *
#* |      To distinguish it with generic 'parInfer' of nl.fitt, arguments changed. | *
#* |         Input:  x, design in nonlinear case is gradient.                      | *
#* |         Input:  theta, without sigma.                                         | *
#* |                 sigma: estimated sigma.                                       | *
#* +-------------------------------------------------------------------------------+ *

pInf <- function(object,confidence=.95){
			.expr1 <- attr(object@predictor,"gradient")							##  F
			.expr2 <- t(.expr1) %*% .expr1										##  F' F
			gtginv <- indifinv( .expr2,F )
			if(is.Fault(gtginv)) gtginv<-ginv(.expr2)
			covmatrix <- object@parameters[["sigma"]] ^ 2 * gtginv			##  sgm2 * (F'F)^-1 
			v2inv <- diag(1/sqrt(diag(covmatrix)))
			corrmat <- v2inv %*% covmatrix %*% v2inv
			parstdev = sqrt(diag(covmatrix ))			

			n <- nrow(.expr1)
			p <- object$form$p
			.expr4 <- sqrt((p+1)*qf(confidence,p+1,n-p-1))
			.expr5 <- parstdev * .expr4
			cilow <- unlist(object$parameters[names(object$form$par)]) - .expr5
			ciupp <- unlist(object$parameters[names(object$form$par)]) + .expr5


			result <- list(covmat=covmatrix ,corrmat = corrmat,parstdev=parstdev,CI=cbind(cilow,ciupp))
			return(result)
	}
