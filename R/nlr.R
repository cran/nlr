#********************************************************************************************************
#  +-------------------------------------------------------------------------------------------------+  #
#**|   Function 'nl' estimating nonlinear regression using all methods.                              |**#
#**|                                                                                                 |**#
#**|  Note: becarefull to using this function when there is not outlier, it                          |**#
#**|    may not work witout outlier, in this case better to use clasic                               |**#
#**|    argumnts:                                                                                    |**#
#**|      formula:      'nl.form', formula, function or expression.                                  |**#
#**|      data:         data, contains dependents and independents,                                  |**#
#**|                    data.frame, list or named matrix it can be.                                  |**#
#**|      start:        starting values, it must contains selstart                                   |**#
#**|                    for nl.form is created.                                                      |**#
#**|      varianceform: variance function, it must be nl.form of variance models                     |**#
#**|      tau:          starting value of tau. if is null the stored value in                        |**#
#**|                    vardnc object of nl.form will be stored.                                     |**#
#**|      control:      is to nlr.control, when error happens and need control                       |**#
#**|                    manually, program try to control errors but other                            |**#
#**|                    un predicted error like log(0) may happens.                                  |**#
#**|         algorithm:   algorithm for computing                                                    |**#
#**|         method:      method of fitting, either official algortitm or control but not both.      |**#
#**|                          but if in method cllasic methods like CME,CLsME, NLLS mention,         |**#
#**|                          automaticale nonrobust will be applied. It get these values:           |**#
#**|                  "default"(is robust MM),"RME","CME","CLSME","RGME","WME","MLE","OLS","TS","RTS"|**#
#**|      derivfree:    T,F use derivative free, default gradient base.                              |**#
#**|      weights:      vector of weights, for weighted regression, if weights present then variance |**#
#**|                    form and correlation form will be ignored.                                   |**#
#**|                                                                                                 |**#
#**|      robustobj:    "nl.form" of robust loss function, or should be presented the name robustform|**#
#**|                                                                                                 |**#
#**|      robustform:   is "nl.form" of a character of name of robust function in "nl.robfuncs"      |**#
#**|                    defaul of these two is "hampel", only when cllasic method mention            |**#
#**|                    will be ignore.                                                              |**#
#**|                                                                                                 |**#
#**|      varianceform: "nl.form" of variance function.                                              |**#
#**|                                                                                                 |**#
#**|      correlation:  autocorrelated error, form of corStruct but not in nlme.                     |**#
#**|      ...            is extra arguments to any of "nl.forms", such as tuning constants to        |**#
#**|                     robustobj, or variancemodel or formula, but must be carefull.               |**#
#**|                                                                                                 |**#
#**|   Important Note: variance is a product function in sigma, i.e.                                 |**#
#**|       varfunc = sigma^2 * h(f,tau), but from variance model "h" must be presented.              |**#
#**|                                                                                                 |**#
#  +-------------------------------------------------------------------------------------------------+  #
#********************************************************************************************************
nlr <- function(formula, data = parent.frame(), start=getInitial(formula,data),  
                control = nlr.control(minlanda=1 / 2 ^ 10, maxiter=25 * length(start)),
                weights=NULL,
                robustobj = NULL,
                                     # *** values in "nl.robfuncs" variable, default hampel ***
                robustform = c("hampel","huber","bisquare","andrew","halph huber","hampel 2","least square"),
                varianceform=NULL,tau=NULL,
                correlation=NULL,covariance=NULL,...) {
  robustform <- match.arg(robustform)
  if(is.null(robustobj)){
    robustobj <- nlr::nl.robfuncs[[robustform]]
  }
  else 
    if(! class(robustobj)=="nl.form") return(Fault(FN=25))
  
  algorithm <- control$algorithm
  method <- control$method
 
  # +-------------------------------------+
  # |   create start                      |
  # |   must be list with parameter names |
  # +-------------------------------------+
# COPIEND FROM NLSQR

  switch(mode(start), 
         "numeric"={.parameters <- names(start)
                    start <- as.list(start)
                    names(start) <- .parameters; },
         "list"={.parameters <- names(start); },
        "NULL"={
          datalist <- as.list(data)
           .parameters <- parameter.names(formula, datalist)
        },
         stop("\"start\" should be numeric or list")
  )
#  if (!length(.parameters)) stop("names for parameters needed, from formula or from start")
  
  #COPIED FROM NLS
#  pnames <- if (missing(start)) {
#    if (!is.null(attr(data, "parameters"))) {
#      names(attr(data, "parameters"))
#    }
#    else {
#      cll <- formula[[length(formula)]]
#      func <- get(as.character(cll[[1L]]))
#      print(func)
#      if (!is.null(pn <- attr(func, "pnames"))) 
#        as.character(as.list(match.call(func, call = cll))[-1L][pn])
#    }
#  }
#  else names(start)
  
  # +-------------------------------------+
  # |   create data                      .|
  # |   must be list with parameter names |
  # +-------------------------------------+
  
  
  # +-----------------------------------+
  # |   create nl.form from expression..|
  # |   when formula is not nl.form     |
  
  switch(class(formula),
      "formula"={
          formula <- convexpr2nlform(formula,names(data),start)
      },
      "function"={
        #formula <- new("nl.form",formula=formula,p=length(start),
        #               inv=NULL,name="User defined",par=start,arguments=list(...),dependent=dependent,
        #               independent=independent,origin=NULL,selfStart=NULL) 
                                                          # ...... latter correct selfstart .......
      },
      "nl.form"={},
      {
        return(Fault(FN=24))  
      }
    )
  
  # +-----------------------------------+
  # | if a weight is given manually then|
  # |    other cases will be ignored.   |
  # +-----------------------------------+

  if(! is.null(weights)) mcase <- 1                    #*** generalized form only weight, other cases ignored
  else
  {
    if(is.null(varianceform) && is.null(correlation)) mcase <- 2      #*** nlmest.NLM, nlsqr
    if(! is.null(varianceform) && is.null(correlation)) mcase <- 3    #*** heterogeneous case
    if(is.null(varianceform) && ! is.null(correlation)) mcase <- 4    #*** autocorrelated case
    if(! is.null(covariance)) mcase <- 5
  }
  switch(control$initials,
    "quantile"={
      #library(quantreg)
      .temp <- nlrq(formula$origin, data=data,start=start, tau=0.5,trace=control$trace)
      start <-as.list(coef(.temp))
    },
    "ols" ={
      .temp <- nlsqr(formula=formula, data=data,start=start,control=control,...)
      start <- .temp$parameters[1:formula$p]
    },
    "lms"={
      fitga <- nl.lmsGA(formula,data,start)
      start <- fitga$parameters
      start[["sigma"]] <- fitga$scale
    }
  )
  switch(as.character(mcase),
         #******************************************
         #***   generalized form only weight  ******
         #******************************************
         "1" = {
           vmat <- diag(weights^2)
           rmat <- diag(1.0/weights)
           
           switch(method,
                 "OLS" = {               ###    clasic, OLS ********
                   if(control$derivfree || algorithm =="Nelder-Mead")
                     result <- nlsnm(formula=formula,data=data,start=start,control=control,vm=vmat,rm=rmat,...)
                   else
                     result <- nlsqr.gn(formula=formula, data=data, start=start,control=control,
                                      vm=vmat,rm=rmat,...)
                   },
                  # "MLE" ={} # for feature version
                  {                      ###   default if in control method="default"  robust NLM or NM  *****
                   if(control$derivfree || algorithm =="Nelder-Mead")
                     result <- nlmest.NM(formula=formula, data=data, start=start,
                                         control=control,robfunc=robustobj,vm=vmat,rm=rmat,...)
                   else{
                      result <- nlmest.NLM(formula=formula, data=data, start=start,
                                          control=control,vm=vmat,rm=rmat,robfunc=robustobj,...)
                      plot(result)
                      if(is.Fault(result))
                            result <- nlmest.WF(formula=formula, data=data, start=start,
                                           control=control,vm=vmat,rm=rmat,robfunc=robustobj,...)
                   }
                  }
                 )
         },
         #******************************************
         #***     nlmest.NLM or OLS           ******
         #******************************************
         "2" = {
           switch(method,
                 "OLS" = {
                   if(control$derivfree || algorithm =="Nelder-Mead")
                     result <- nlsnm(formula=formula,data=data,start=start,control=control,...)
                   else{
                     
                     result <- nlsqr(formula=formula, data=data,start=start,control=control,...)
                   }
                 },
                 "lms"={
                   fitga <- nl.lmsGA(formula,data,start)
                   result <- fitga
                 },
                 {          #**** default: is equivalent to NLM   ****
                   if(control$derivfree || algorithm =="Nelder-Mead")              #***** derivative free ********
                      result <- nlmest.NM(formula=formula, data=data, start=start,
                                        control=control,robfunc=robustobj,...)
                   else {                               #*****  gradient base  ********
                     if(control$singularCase == 1) {
                        result <- nlmest.NLM(formula=formula, data=data, start=start,
                                         robfunc=robustobj,control=control)
                        if(is.Fault(result)){
                          result <- nlmest.WF(formula=formula, data=data, start=start,
                                           robfunc=robustobj,control=control)
                        }
                        if(is.Fault(result)){
                           result <- nlmest.NLM.sCase2(formula=formula, data=data, start=start,
                                                   robfunc=robustobj,control=control)
                        }
                     }
                     else{
                       result <- nlmest.NLM.sCase2(formula=formula, data=data, start=start,
                                            robfunc=robustobj,control=control)
                       if(is.Fault(result))
                          result <- nlmest.WF(formula=formula, data=data, start=start,
                                             robfunc=robustobj,control=control)
                     }
                   }
                 }
           )
         },
         #******************************************
         #***        heterogeneous case        *****
         #******************************************
         "3" = {
           
           if(control$derivfree || algorithm =="Nelder-Mead")              #***** derivative free ********
           {
             switch(method,
                    "RME" = {               #***** default from nlr.control ***
                      result <- dfr.robhetro(formula=formula,data=data,start=start,tau=tau,
                                            robfunc=robustobj,varmodel=varianceform,control=control,method="NM",...)
                    },
                    "RGME" = {
                      result <- dfr.robhetroLS(formula=formula,data=data,start=start,robfunc=robustobj,
                                              tau=tau,varmodel=varianceform,...)
                    },
                    "WME" = {
                      ### must be developed yet, WME is deriv free
                      result <- nl.robhetroWM(formula=formula,data=data,start=start,robfunc=robustobj,
                                              tau=tau,varmodel=varianceform,control=control,...)
                    },
                    "CME" = {
                      result<- dfr.hetro(formula=formula, data=data, start=start,
                                        tau=tau,varmodel=varianceform,...)
                    },
                    "CLSME" = {
                      result <- dfr.hetroLS(formula=formula, data=data, start=start,
                                           tau=tau,varmodel=varianceform,...)
                    },
                    {          #**** default: is equivalent to RME above   ****
                      result <- dfr.robhetro(formula=formula,data=data,start=start,tau=tau,
                                             robfunc=robustobj,varmodel=varianceform,control=control,...)  
                    }
             )
           }
           else                               #*****  gradient base  ********
           {
             
              switch(method,
                     "RME" = {               #***** default from nlr.control ***
                        result <- nl.robhetro(formula=formula,data=data,start=start,tau=tau,
                                              robfunc=robustobj,varmodel=varianceform,control=control,...)
                      },
                     "RGME" = {
                        result <- nl.robhetroLS(formula=formula,data=data,start=start,robfunc=robustobj,
                                              tau=tau,varmodel=varianceform,...)
                      },
                     ### must be developed yet, WME is deriv free
                     "WME" = {
                        result <- nl.robhetroWM(formula=formula,data=data,start=start,robfunc=robustobj,
                                             tau=tau,varmodel=varianceform,control=control,...)
                      },
                     "CME" = {
                       
                       result<- nl.hetro(formula=formula, data=data, start=start,
                                      tau=tau,varmodel=varianceform,...)
                     },
                     "CLSME" = {
                       result <- nl.hetroLS(formula=formula, data=data, start=start,
                                        tau=tau,varmodel=varianceform,...)
                     },
                     {            # ***** default: is same as RME above ****
                       result <- nl.robhetro(formula=formula,data=data,start=start,tau=tau,
                                             robfunc=robustobj,varmodel=varianceform,control=control,...)
                     }
              )
           }
         },
         #******************************************
         #***     autocorrelated case          *****
         #******************************************
        "4" = {
          switch(method,
                 "TS" =,"OLS"= {
                   if(control$derivfree || algorithm =="Nelder-Mead"){
                     result <-  dfr.corrts(formula=formula,data=data,start=start,control=control,correlation=correlation,...)
                   }
                   else{
                     result <-  nl.corrts(formula=formula,data=data,start=start,control=control,correlation=correlation,...)
                   }
                 },
                 "RTS" =, "default" = {
                   if(control$derivfree || algorithm =="Nelder-Mead"){
                     resutlt <- dfr.robcorrts(formula=formula, data=data, start=start, 
                                               control=control,correlation=correlation,robfunc=robustobj,...)
                   }
                   else{
                     result <- nl.robcorrts(formula=formula, data=data, start=start, 
                                             control=control,correlation=correlation,robfunc=robustobj,...)
                   }
                 }
                 )
         },
         #******************************************
         #***     Generalized, manual covariance ***
         #******************************************
         "5" = {
           vmat <- covariance
           rmat <- chol(vmat)
           switch(method,
                  "OLS" = {               ###    clasic, OLS ********
                    if(control$derivfree || algorithm =="Nelder-Mead")
                      result <- nlsnm(formula=formula,data=data,start=start,control=control,vm=vmat,rm=rmat,...)
                    else
                      result <- nlsqr.gn(formula=formula, data=data, start=start,control=control,
                                         vm=vmat,rm=rmat,...)
                  },
                  {                      ###   default if in control method="RME"  robust NLM or NM  *****
                    if(control$derivfree || algorithm =="Nelder-Mead")
                      result <- nlmest.NM(formula=formula, data=data, start=start,
                                          control=control,robfunc=robustobj,vm=vmat,rm=rmat,...)
                    else
                      result <- nlmest.NLM(formula=formula, data=data, start=start,
                                           control=control,vm=vmat,rm=rmat,robfunc=robustobj,...)
                  }
                  )
         }
         )
  if(is.Fault(result)){
    if(class(result)=="Fault"){
      rslt <- nl.fitt(
            form =         formula,
            data =         as.list(data),
            sourcefnc =    match.call(),
            Fault =        result
            )
      result <- rslt
    }
    else{
      rslt <- nl.fitt(
        form =         formula,
        data =         as.list(data),
        sourcefnc =    match.call(),
        Fault =        result@Fault
      )
      result <- rslt
    }
      
  }
  result@sourcefnc = match.call()
  return(result)
           
}
#+#################################################################################+
#|                                                                                 |
#|                         End of the nlr                                          |
#|                                                                                 |
#|                     developed for R 23/11/2013                                  |
#|                                                                                 |
#|                     Hossein Riazoshams, Department of Statistics                |
#|                                                                                 |
#|                          Stockholm university                                   |
#+#################################################################################+

##### example #####
#b5=list(xr=Weights$Date,yr=Weights$Weight)
#start<-list(p1=2200,p2=38,p3=.1)

#aa1 <- nlr(formula=nlrobj1[[14]],data=b5,start=start,robustobj=nl.robfuncs[["hampel"]],
#                     tau=list(sg=.09,landa=2),varianceform=nlrobjvarmdls1[[1]],control=nlr.control(tolerance=1e-8,
#                   maxiter=500,derivfree=F),method="CME")#,delta=c(0.2,1,1,160,.2,1,.03))


#nlr(y~a+exp(-b*x),start=list(a=1,b=3),data=list(x=c(1,2,3),y=c(2,3,4)),
#    control=nlr.control(method="ML"),robustform="hamp")  
  
#method=c("RME","CME","CLSME","RGME","WME","MLE","OLS"))



