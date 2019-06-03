library(MASS)
# +-----------------------------------------------------------------+ #
# |  convexpr2nlform create nl.form from an expression              | #
# |    arguments:                                                   | #
# |       form: two sided or one sided or function.                 | #
# |       namesdata: name of data include independent and possibly  | #
# |                  dependent in two sided fomula.                 | #
# |       start: is list of parameters.                             | #
# |   result is nl.form object with gradient and hessian also.      | #
# +-----------------------------------------------------------------+ #
convexpr2nlform<-function(form,namesdata=NULL,start,inv=NULL,name="User Defined",...){
  n <- length(form)
  if(n==2){
    a0 <- call("~", form[[2]])
    a1 <- formula(a0)
    a2 <- names(start)
    b <- deriv3(a1, a2)
    depname <- NULL
    formula <- call("~", b)
    frmnames <- all.vars(form[[2]])
    pnames <- names(start)
    matchnames<-charmatch(pnames,frmnames)
    indepnames <- frmnames[-matchnames]
  }
  if(n==3){
    a0 <- call("~", form[[3]])
    a1 <- formula(a0)
    a2 <- names(start)
    b <- deriv3(a1, a2)
    depname <- all.vars(form[[2]])     # name of dependent variable
    d2 <- call(depname)
    formula <- call("~",d2[[1]], b)
    frmnames <- all.vars(form[[3]])
    pnames <- names(start)
    matchnames<-charmatch(pnames,frmnames)
    indepnames <- frmnames[-matchnames]
  }
  new("nl.form",formula=formula,formtype="formula",p=length(start),
      name=name,par=start,dependent=depname,
      independent=indepnames,origin=form,inv=inv,...)
                              ## corect self start later ##
}

# +-----------------------------------------------------------------+ #
# |                       Hossein Riazoshams                        | #
# |                          2014/04/19                             | #
# |            Department of statistics, stockholm university       | #
# +-----------------------------------------------------------------+ #
