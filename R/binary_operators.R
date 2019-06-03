#*********************************************************************************
#**   function: %c%, implementing crossproduct in splus                         **
#**      compute: t(x) * y=crossprod(x,y)                                       **
#**      variables:  x (m*n), y(m*q)                                            **
#**      result:  is (n*q) matrix    =   sum(x[,i]*y[,j]).                       **
#*********************************************************************************

"%c%" <- function(x,y){ crossprod(x,y)}

#*********************************************************************************
#**   function: 3d*m, product array to matrix transform                         **
#**      compute: 'array * vector', which is (n,p,p) array.                     **
#**               matrix multiple each vector of arry[,i,j] is n.1              **
#**      variables:  ary: is (n*p*p) array.                                     **
#**      				vector: is a matrix with dimension (n*1)                        **
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
#**      result:  is (n*p*p) hessian,                                           **
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
