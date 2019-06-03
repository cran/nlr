#*********************************************************************************
#**   M-estimate of scale                                                       ** 
#*********************************************************************************
mscale<-
function(u)
{
	# Scale M-estimator with 50% breakdown
	# Yohai (1987) Annals, Stromberg (1993) JASA.
	#
	# GKS  2 June 99
	#
	if(mean(u == 0) >= 0.5) return(0)
	U <- abs(u)
	s <- median(U)/0.6744898
	iter <- 0
	repeat {
		iter <- iter + 1
		z <- u/0.212/s
		d1 <- mean(rho.hampel(z)) - 3.75
		d2 <- mean(z * psi.hampel(z))
		s <- s * (1 + d1/d2)
		if(iter > 50) {
			cat("mscale: Max iterations exceeded")
			break
		}
		if(abs(d1/d2) < 1e-014)
			break
	}
	s
}
