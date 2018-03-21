#*********************************************************************************
#**   function: prodVAV product array to matrix                    				**
#**      compute: 'vector (T) * array * vector', which is (n,1) vector.			**
#**      variables:  ary: is (p*p*n) array.										**
#**		  				vector: is a vector with dimension (p)					**
#*********************************************************************************
prodVAV <- function(ary, vector)
{
	p <- length(ary[, 1, 1])
	n <- length(ary[1, 1, ])
	prd <- array(rep(0, times=n), dim=n)
	for (k in 1:n)
	{
		t1 <- ary[, , k]
		t2 <- vector %*% t1
		t4 <- t(vector)
		t3 <- t2 %*% vector
		prd[k] <- t3
	}
	
	return(prd)
}

#*********************************************************************************
#**   function: prodAV product array to matrix                    				**
#**      compute: 'array * vector', which is (n,p) vector.						**
#**      variables:  ary: is (p*p*n) array.										**
#**		  				vector: is a vector with dimension (p)					**
#*********************************************************************************
prodAV <- function(ary, vector)
{
	p <- length(ary[, 1, 1])
	n <- length(ary[1, 1, ])
	prd <- array(rep(0, times=n * p), dim=c(n, p))
	for (k in 1:n)
	{
		t1 <- ary[, , k]
		t2 <- t1 %*% vector
		prd[k, ] <- t2
	}
	
	return(prd)
}
#*********************************************************************************
#**   function: 3d*v, product array to matrix                    				**
#**      compute: 'array * vector', which is (n,p,p) vector.					**
#**      variables:  ary: is (n*p*p) array.										**
#**		  				vector: is a vector with dimension (n)					**
#*********************************************************************************
"%3d*v%" <- function(ary, vector)
{
	p <- length(ary[ 1, 1,])
	n <- length(ary[,1, 1 ])
	prd <- array(0,c(n,p,p))
	for (k1 in 1:p)
		for(k2 in 1:p)
		{
			t1 <- ary[,k1 ,k2 ]
			t2 <-  vector * t1 
			prd[,k1,k2] <- t2
		}
	return(prd)

}
#*********************************************************************************
#**   function: 3d*m, product array to matrix transform                         **
#**      compute: 'matrix * vector', which is (n,p,p) array.                    **
#**               matrix multiple each vector of arry[,i,j] is n.1              **
#**      variables:  ary: is (n*p*p) array.                                     **
#**		  				vector: is a matrix with dimension (n*n)                    **
#*********************************************************************************
"%3d*m%" <- function(ary, vector)
{
	p <- length(ary[ 1, 1,])
	n <- length(ary[,1, 1 ])
	prd <- array(0,c(n,p,p))
	for (k1 in 1:p)
		for(k2 in 1:p)
		{
			t1 <- ary[,k1 ,k2 ]
			t2 <-  sum(vector * t1)
			prd[,k1,k2] <- t2
		}
	return(prd)

}

#*********************************************************************************
#**   function: m3d, matrix 3 dimentional product.                              **
#**      compute: 'm (p*n) * (n*p)'                                             **
#**      variables:  mat1 (p*n), mat2(n*p), (both like gradient)                **
#**		  result:  is (n*p*p) hessian,                                           **
#**          which is n: mat1[,i] %*% mat2[i,], (p.1) * (1.p) = (p.p), i=1..n   **
#*********************************************************************************

"%m3d%" <- function(mat1, mat2)
{
	n1 <- dim(mat1)[2]
	p1 <- dim(mat1)[1]
	n2 <- dim(mat2)[1]
	p2 <- dim(mat2)[2]
	if((n1!=n2) || (p1!=p2)) stop("Problem in '%m3d%.default'(b, a): Number of rows of x should be the same as number of 
columns of y, and Number of columns of x should be the same as number of rows of y, Use traceback() to see the call stack")

	prd <- array(0,c(n1,p1,p1))
	for(i in 1:n1){
		m1 <- mat1[,i]
		m2 <- mat2[i,]
		prd[i,,] <- m1 %*% t(m2)
	}
	return(prd)

}

#*********************************************************************************
#**   function: prodAV product array to matrix                                  ** 
#**      compute: 'array * vector', which is (p,p) matrix.                      **
#**      variables:  ary: is (n*p*p) array. Diffrent from before, then change.  **
#**		  				vector: is a vector with dimension (n)                     **
#*********************************************************************************
prodVA <- function(ary, vector)
{
	if(class(ary)=="numeric") return(sum(ary*vector))
	p <- length(ary[ 1, 1,])
	n <- length(ary[,1, 1 ])
	prd <- matrix(rep(0, times=p * p),nrow=p)
	for (k1 in 1:p)
		for(k2 in 1:p)
		{
			t1 <- ary[,k1 ,k2 ]
			t2 <-  sum(vector * t1) 
			prd[k1,k2] <- t2
		}
	return(prd)
}

#------------------------------------------------------------------------------------------
#*****   compute QR
#------------------------------------------------------------------------------------------
fullqr <- function(x)
{
	qrl <- qr(x)
	r <- qrl$qr
	r[row(r) > col(r)] <- 0
	
	# make lower triangular
	u <- qrl$qr
	u[row(u) < col(u)] <- 0
	
	# make upper triangular
	u[row(u) == col(u)] <- qrl$qraux
	
	# replace diagonal elements
	q <- diag(nrow(u))
	for (j in 1:ncol(u))
	{
		h <- diag(nrow(u)) - outer(u[, j], u[, j]) / qrl$qraux[j]
		q <- q %*% h
	}
	
	n <- nrow(q)
	p <- qrl$rank
	q1 <- q[c(1:n), c(1:p)]
	
	# partition Q=[q1|q2]
	q2 <- q[c(1:n), (p + 1):n]
	
	# partition Q=[q1|q2]
	r1 <- r[c(1:p), c(1:p)]
	
	# partition R=[R1/0]
	rinv <- ginv(r1)
	result <- list(q=q, r=r, q2=q2, r1=r1, rinv=rinv, qr=qrl)
	return(result)
}

#*****************************************************************************************************
individ <- function(hessian)
{
	n <- length(hessian[, 1, 1])
	p <- length(hessian[1, 1, ])
	sp <- matrix(hessian[, 1, 1], ncol=1)
	for (i in 2:p)
	{
		sp <- cbind(sp, hessian[, i, 1:i])
	}
	return(sp)
}


#*********************************************************************************
#**   function: %c%, implementing crossproduct in splus                         **
#**      compute: t(x) * y=crossprod(x,y)                                       **
#**      variables:  x (m*n), y(m*q)                                            **
#**  	  result:  is (n*h) matrix    =   sum(x[,i]*y[,j]).                       **
#*********************************************************************************

"%c%" <- function(x,y) crossprod(x,y)