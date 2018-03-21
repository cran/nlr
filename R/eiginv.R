#*****************************************************************************
#* +-----------------------------------------------------------------------+ *
#* |   eiginv: function to compte inverse of a matrix using eigen values   | *
#* |           and eigenvectors.                                           | *
#* |   mtrx:       matrix for computing.                                   | *
#* |   symmetric:  better to be entred by user as T if the matrix is       | *
#* |               symetric due to computer memmory significan digit.      | *
#* |   stp:        if true the procedure stop in singular case.            | *
#* |  THe inverse of a matrix can be expressed by eigenvalues and eigenve_ | *
#* |  ctors of a matrix as:                                                | *
#* |      A^-1 = sum of[ (1/landa i) * ei * ei(T) ]                        | *
#* |           = E * diag(1/landa'i) * E T , as matrix form is faster.     | *
#* |                                                                       | *
#* |  A symmetric matrix (A) is psitive definit if and only if every       | *
#* |     eigenvalues of (A) s positive.                                    | *
#* |  A is nonnegatie definite matrix if and only if all of its eigenvalues| *
#* |     are greather than or equal to zero. see Johnson and Whichern      | *
#* |  page 63.                                                             | *
#* |                                                                       | *
#* |    After 14 years of my bachelore, my programing brain reactivated    | *
#* |                                21-oct-2008                            | *
#* |                             Hossein Riazoshams                        | *
#* +-----------------------------------------------------------------------+ *
#*****************************************************************************

eiginv <- function(mtrx,stp=T,symmetric=all(mtrx==t(mtrx))){
	if(any(is.na(mtrx))|| any(is.inf(mtrx))) return(Fault(FN=12))
	n <- nrow(mtrx)
	if(n==1) return(1.0/mtrx)
	mtrxeig <- eigen(mtrx,symmetric=symmetric)
	if( ! ( all( mtrxeig$values > 0) ) & symmetric ){
		if(stp){
			stop("Matrix is not possitive definite, in eiginv")
		}
		return(Fault(FN=9))
	}
	rvs <- 1.0 / mtrxeig$values
	.expr1 <- diag(rvs)
	.expr2 <- .expr1 %*% t(mtrxeig$vectors)
	mtrxinv <- mtrxeig$vectors %*% .expr2
	return(mtrxinv)
}
###############################
indifinv <- function(mtrx,stp=T,symmetric=all(mtrx==t(mtrx))){
	if(any(is.na(mtrx)) || any(is.inf(mtrx))) return(Fault(FN=12,FT="indifinv"))
  n <- nrow(mtrx)
  if(n==1) return(1.0/mtrx)
  mtrxeig <- eigen(mtrx,symmetric=symmetric)
  rvs <- 1.0 / mtrxeig$values
	if(any(is.inf(rvs)) ){ 
		if(stp){
			stop("Matrix is not possitive definite, in indifinv")
		}
		return(Fault(FN=9,FF = "indifinv"))
	}
  
	.expr1 <- diag(rvs)
	.expr2 <- .expr1 %*% t(mtrxeig$vectors)
	mtrxinv <- mtrxeig$vectors %*% .expr2
  return(mtrxinv)
}



