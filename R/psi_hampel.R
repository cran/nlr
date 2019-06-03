
psi.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  2 June 99
#
	U <- abs(u)
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	psi <- u
	psi[B] <- sign(u[B]) * a
	psi[C] <- sign(u[C]) * a * (c - U[C])/(c - b)
	psi[D] <- 0
	psi
}

