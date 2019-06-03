#**********************************************************************************************
#**   Robust Two stage estimate (RTS)                                                        **
#**  Arguments:                                                                              **
#**   correlation: correlation structures not follow nlme package, but argument pass to nlme. **
#**     after giving StructName the argument to appropriate nlme must be given               **        
# **     list of arguments:                                                                  **
#**        StructName: character name of structure, which pass to nlme, corAR1, corARMA,     **
#**                                         corCAR1, corCompSymm, corExp, corGaus, corLin,   **
#**                                         corRatio, corSpher, corSymm                      **
#**       if correlation$StructName = above the following argument to TS function must provide**
#**        "corCAR1",                                                                        **
#**                   pass to  "corAR1(value, form, fixed)", value=p compute from            **
#**                            robar(ri, order = 1), other argument see ?corAR1              **
#**        "corARMA",                                                                        **
#**                   pass to  "corARMA(value, form, p, q, fixed)", value=p compute from     **
#**                            tm=arimax, tm$coef pass to "value" for phi& theta, as p & q   **
#**                            p,q comes from correlation(order(p,d,q), seasonal not set yet)**
#**       "manual" manualcorrelation, a function for computing autocorrelation must be       **
#**                   given by user in correlation$manualcorr,                               **
#**                   which is name as character of direct function                          **
#**       others not defined yet.                                                            **
#**       order and seasonal, and other arguments are in ... pass to robust arimax function  **
#** temporarty correlation parameters like order, seasonla and so on goes into ... argument  **
#**********************************************************************************************
nl.robcorrts <-
  function(formula,
           data,
           start = getInitial(formula, data),
           control = nlr.control(
             tolerance = 0.0010,
             minlanda = 1 / 2 ^ 10,
             maxiter = 25 * length(start)
           ),
           correlation = list(
             StructName = "NAN",
             manualcorr = NULL
           ), 
           robfunc,
           ...)
  {
    q <- correlation
    #	tols2 <- nl.mm3th(formula, data=data, start=start, nl.robfuncs[[2]])
    tols2 <-
      nlmest.NLM(
        formula,
        data = data,
        start = start,
        robfunc = robfunc,
        control = control,
        ...
      )
    if (is.Fault(tols2))
      return(tols2)
    ri <- residuals(tols2)
    n <- length(ri)
    dotlist <- list(...)
    switch(
      correlation$StructName,
      "corAR1" = {
        tm <- robar(ri, order = 1)
        cs1<-corAR1(tm[[2]],form=~1)
        cs2 <- Initialize(cs1, data = as.matrix(ri))
        vmat <- corMatrix(cs2)
        #vinv <- solve(vmat)
        #umat <- chol(vmat)
        #ut <- t(umat)
        #rmat <- solve(ut)
        
        .temp1 <- 1.0 / (1.0 - tm$ar ^ 2)
        vinv <- diag(c(1, rep(1 + tm$ar ^ 2, n)))
        vinv[col(vinv) == row(vinv) + 1] <- -tm$ar
        vinv[row(vinv) == col(vinv) + 1] <- -tm$ar
        rmat <- diag(c(sqrt(1 - tm$ar ^ 2), rep(1, n - 1)))
        rmat[row(rmat) == col(rmat) + 1] <- -tm$ar
        rmat <- rmat / sqrt(1 - tm$ar ^ 2)
        autpar <- tm$ar
      },
      
      "corARMA" = {
        arimalist<- as.list(formals(arimax))
        matcharimax <- match.arg(names(correlation),names(arimalist),several.ok = TRUE)
        arimaxarguments<-correlation[matcharimax]
        pcorr <- correlation$order[1]
        qcorr <- correlation$order[3]
        ncorr <- pcorr + qcorr
        #tm <-arimax(ri,order = order,include.mean = FALSE)
        tm <- do.call(arimax,c(list(x=ri),arimaxarguments))
        correst <-corARMA(tm$coef[1:(pcorr+qcorr)],p = pcorr,q = qcorr,form=~1,fixed=~1)
        cs2 <- Initialize(correst, data = as.matrix(ri))
        vmat <- corMatrix(cs2)
        v2 <- eiginv(vmat, symmetric = T, stp = F)
        if (is.Fault(v2))
          return(v2)
        for (i in 1:n)
          for (j in i:n)
            v2[i, j] <- v2[j, i]
        rmat <- chol(v2)
        autpar <- tm$coef
      },
      "corCAR1" = {
        
      },
      "corCompSymm" = {
        
      },
      
      "corExp" = {
        
      },
      "corGaus" = {
        
      },
      "corLin" = {
        
      },
      "corRatio" = {
        
      },
      "corSpher" = {
        
      },
      "corSymm" = {
        
      },
      "manual" = {
        # .... manual function must have ri as entry 
        #..... and return n*n matrix of autocorrelation function
        ftar<-do.call(correlation$manualcorr,list(ri=ri))
        tm<-ftar$tm
        vmat<-ftar$vmat
        v2 <- eiginv(vmat, symmetric = T, stp = F)
        if (is.Fault(v2))
          return(v2)
        for (i in 1:n)
          for (j in i:n)
            v2[i, j] <- v2[j, i]
        rmat <- chol(v2)
        autpar <- tm$coef
      },
      "NAN" = {
        vmat<-NULL
        rmat<-NULL
      }
      
    )
    tolerance <- control$tolerance * 1e3
    minlanda <- control$minlanda / 1e4
    t2st <- nlmest.NLM(
      formula,
      data = data,
      start = tols2$parameters[names(formula$par)],
      # withought sigma
      robfunc = robfunc,
      vm = vmat,
      rm = rmat,
      control = nlr.control(tolerance = tolerance, minlanda =
                              minlanda),
      ...
    )
    if (is.Fault(t2st))
      return(t2st)
    t2st@autpar <- as.list(autpar)
    t2st@autcorr <- list(tm = tm)
    t2st@method@subroutine <- "nl.robcorrts"
    t2st@method@subroutineBR <- "nlmest.NLM"
    t2st@method@method <- "RTS"
    return(t2st)
  }