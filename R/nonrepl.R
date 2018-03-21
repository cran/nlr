
#****** data$x and find measures for data$y
nonrepl <- 
function(data)
{
	# marks replications only if there are more than 2
	# (otherwise use of z_values is preferable to use of s_i`s
	# to find replications in the design
	x <- data$x
	if(!is.null(dim(x)))
		stop("This version allows only\none exploratory variable")
	o <- order(x)
	x <- x[o]
	n <- length(x)
	z <- x[diff(x) > 0]
	k <- length(z)
	if(z[k] < x[n])
		z <- c(z, x[n])
	xm <- match(z, x)
	xk <- x[xm]
	ni <- c(diff(xm), n - xm[length(xk)] + 1)
	data$x <- x
	data$y <- data$y[o]
	if(!is.null(data$fw))
		data$fw <- data$fw[o]
	if(!is.null(data$res))
		data$res <- data$res[o]
	k <- length(z)
	# now start to extend points with ni==2
	k2 <- k + sum(ni == 2)
	ni2 <- numeric(k2)
	xk2 <- numeric(k2)
	xm2 <- numeric(k2)
	xo <- numeric(n)
	j <- 1
	l <- 0
	for(i in 1:k) {
		if(ni[i] == 2) {
			ni2[j] <- ni2[j + 1] <- 1
			xk2[j] <- xk2[j + 1] <- xk[i]
			xm2[j] <- xm[i]
			xm2[j + 1] <- xm[i] + 1
			xo[l + 1] <- j
			xo[l + 2] <- j + 1
			j <- j + 2
		}
		else {
			ni2[j] <- ni[i]
			xm2[j] <- xm[i]
			xk2[j] <- xk[i]
			for(m in 1:ni[i])
				xo[l + m] <- j
			j <- j + 1
		}
		l <- l + ni[i]
	}
	data$xk <- xk2
	data$ni <- ni2
	data$xm <- xm2
	data$k <- k2
	data$xo <- xo
	yq <- numeric(k2)
	ys <- numeric(k2)
	for(i in 1:k2){
		yq[i] <- mean(data$y[xo == i])
		ys[i] <- sd(data$y[xo == i])
	}
	data$yq <- yq
	data$ys <- ys
	data
}
