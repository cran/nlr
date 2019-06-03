
#*********************************************************************************
#**   function: prodVA product array to matrix                                  ** 
#**      compute: 'array * vector', which is (p,p) matrix.                      **
#**      variables:  ary: is (n*p*p) array. Diffrent from before, then change.  **
#**  	  				vector: is a vector with dimension (n)                     **
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
