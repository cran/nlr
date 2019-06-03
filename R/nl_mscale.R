#************************************************************************
#**   nl.mscale function                                               **

#**	# Scale M-estimator with 50% breakdown                              **
#**	# Yohai (1987) Annals, Stromberg (1993) JASA.                       **
#**	#                                                                   **
#**	# GKS  2 June 99                                                    **
#**	#                                                                   **
#**     use newton raphson to solve                                    **
#**             (1/n) sum(rho(ri/sgm)) = b                             **
#**     by modified itterated weighted least square.                   **
#**           sgm(n+1) = sgm * [ 1 + d1/d2 ]                           **
#**           d1 = 1/n sum(rho(ri/sgm*k1))                             **
#**           d2 = 1/n sum( ri/k.sgm psi(ri/ksgm) )                    **
#**                                                                    **
#**          Re aranged by Hossein Riazoshams.                         **
#**                 23/09/2009                                         **
#************************************************************************
nl.mscale <- function(u,robfunc,...)
{
	if(mean(u == 0) >= 0.5) return(0)
	U <- abs(u)
	s <- median(U)/0.6744898
	iter <- 0
	repeat {
		iter <- iter + 1
		z <- u/robfunc$arguments$k0/s
		rhof <- robfunc$fnc(z,...)
		d1 <- mean(as.numeric(rhof))  - robfunc$arguments$maxrho5
		d2 <- mean(z * attr(rhof,"gradient"))
		s <- s * (1 + d1/d2)
		if(iter > 50) {
			return(Fault(FN=1,FF="nl.mscale"))
		}
		if(abs(d1/d2) < 1e-014)
			break
	}
	return(s)
}