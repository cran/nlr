rho.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Integral of Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  31 May 99
#
	U <- abs(u)
	A <- (U <= a)	#increasing
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	rho <- U
	rho[A] <- (U[A] * U[A])/2
	rho[B] <- a * (U[B] - a/2)
	rho[C] <- a * (b - a/2) + a * (U[C] - b) * (1 - (U[C] - b)/(c - b)/2)
	rho[D] <- (a * (b - a + c))/2
	rho
}

