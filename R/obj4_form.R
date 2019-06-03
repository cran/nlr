####################################################################################
#**     Set of nonlnear classes, (nl).
#**     formula can be function, and call, expression and formulas are 
#**     known as call
#**     classes: nl.form, nl.fitt

#+#################################################################################+
#**       nl.form:  a nonlinear form, for storing formulas and paameters...         |
#**		Methods:                                                               \
#**			nl.form: constructor function, decide aboute the type of form,   |
#**				  fomula, (expression & call) changes to formula, for    /
#**				  generallity in feature it may have to be change.       |
#**			eval:  generic function to evaluate the form.                  / 
#**	Writhen by Hossein Riazoshas, Date: 24/09/2008, Hour: 2am                /
#**
#**
#+--------------------------------------------------+
#|      class:  nl.form                             |
#+--------------------------------------------------+


setClass("nl.form",representation(
									formula=         "callorNULL",
									fnc=             "functionorNULL",
									dependent =      "characterorNULL",
									independent =    "characterorNULL",
									formtype =       "character",
									p =              "numericorNULL",
									inv =            "callorNULL",
									name =           "character",
									par =            "list",
									arguments =      "list",         ## other parameters like a,b... 
									origin=          "callorNULL",									
									selfStart=       "functionorNULL"
									                                ##  in robust functions
									)
			)

setClassUnion("nl.formorNULL", c("nl.form", "NULL"))
###################################################
##   constructor
###################################################
nl.form<-function(form,p=NULL,inv=NULL,name="nl.frm",par=NULL,dependent=NULL,independent=NULL,origin=NULL,selfStart=NULL,...){
	datalist <- as.list(par)
	if(is.null(dependent)) print("dependent missed")
	switch( mode(form), 
		"call"={
			fm1=(form)
			fm2 = NULL
			ft = "formula";},
		"expression"={
			fm1=as.formula(form)
			fm2 = NULL
			ft = "formula";},
		"formula"={
			fm1=form
			fm2 = NULL
			ft = "formula";},
		"function"={
			fm1 = NULL
			fm2 = form
			ft = "function"; },
    {
      return(Fault(FN=24,FF="nl.form"))
    }
	)

  new("nl.form",formula=fm1,fnc=fm2,formtype=ft,p=p,
		inv=inv,name=name,par=datalist,arguments=list(...),dependent=dependent,
		independent=independent,origin=origin,selfStart=selfStart)
}

#+#################################################################################+
#|            Method (eval) generic. to evaluate the formula at data list.         |
#|           make decisionabout the formula type, call, function formula or        |
#|           expression.                                                           |
#|            result is computed real value                                        |
#|            different from splus, with no ... argument                           |
#+#################################################################################+
setMethod("eval",
    signature(expr = "nl.form"),
    function (expr, envir = parent.frame(), enclos = if (is.list(envir) || 
        is.pairlist(envir)) parent.frame() else baseenv()) 
    {
		.datalist <- as.list(envir)
		result <- NULL
		switch(expr@formtype,
			"formula" = { 
					switch(as.character(length(expr@formula)), 
					"1" = {result$response = NULL
						if(mode(expr@formula[[2]])=="call") result$predictor = eval(as.formula(expr@formula[[2]]), envir)
						else result$predictor = eval(expr@formula[[2]], envir); },

					"2" = {result$response = NULL
						if (mode(expr@formula[[2]]) =="call") result$predictor = evald(as.formula(expr@formula[[2]]), envir)
						else result$predictor = eval(expr@formula[[2]], envir); },

					"3" = { if (mode(expr@formula[[2]]) =="call") result$response = eval(as.formula(expr@formula[[2]]), envir)
							else result$response = eval(expr@formula[[2]], envir)			##   left hand formula
              if (mode(expr@formula[[3]]) =="call") result$predictor = eval(as.formula(expr@formula[[3]]), envir)
							else result$predictor = eval(expr@formula[[3]], envir);
						},			##   right hand formula
					 result=NULL); 
			},
																			###############################################
																			##  In function case, the user must enter the##
																			##  arguments of function into a list with   ##
																			##  names of function argument inside it, not##
																			##  not arguments seperatlly.                ##
			"function"={										##  function case                            ##
				.nf <- length(expr@fnc)       ##  length of all argments and function body.##
				.arglist <- expr@arguments		##  default arguments of the object in slot  ##
																			##                               arguments   ##
				temp <- formals(expr@fnc)   ##  list of argumetns and default values     ##
				.arglist[names(temp)] <- temp ##  A list of arguments and default values,  ##
																			## and object slot 'arguments' if not changed##
				.arglist[names(envir)] <- envir  ##  Substitude data to arguments, new values ##
																			##  will be replaced in default values and   ##
																			##  non default values authomatically replaced#
																			##  by called values.                        ##
				result <- do.call(expr@fnc,.arglist)
			}, 
			result=NULL
		)
 		if(! is.atomic(result)){
      if(! is.null(result$response))                          ##  if any of the response or predictor be   ##
		    if(any(is.nan(result$response))) result=Fault(FN=18)  ##  'nan' , NULL value return back.          ##
		  if(! is.null(result$predictor))                         ##  due to any reason it can be, zero devide ##
  		  if(any(is.nan(result$predictor))) result=Fault(FN=18) ##  missing values ...., in this case feature##
		}
    return(result)															   ##  missing values ...., in this case feature##
    }
)


#+#################################################################################+
#|            Method (evald) generic. to evaluate the formula at data list.        |
#|           make decisionabout the formula type, call, function formula or        |
#|           expression.                                                           |
#|           This is extension to eval to embed ... argument                       |
#|            result is computed real value                                        |
#+#################################################################################+

eval.nl.form=function(expr,envir,...)    {
  .datalist <- as.list(envir)
  result <- NULL
  switch(expr@formtype,
         "formula" = { 
           switch(as.character(length(expr@formula)), 
                  "1" = {result$response = NULL
                         if(mode(expr@formula[[2]])=="call") result$predictor = eval(as.formula(expr@formula[[2]]), envir)
                         else result$predictor = eval(expr@formula[[2]], envir); },
                  
                  "2" = {result$response = NULL
                         if (mode(expr@formula[[2]]) =="call") result$predictor = evald(as.formula(expr@formula[[2]]), envir)
                         else result$predictor = eval(expr@formula[[2]], envir); },
                  
                  "3" = { if (mode(expr@formula[[2]]) =="call") result$response = eval(as.formula(expr@formula[[2]]), envir)
                          else result$response = eval(expr@formula[[2]], envir)			##   left hand formula
                          if (mode(expr@formula[[3]]) =="call") result$predictor = eval(as.formula(expr@formula[[3]]), envir)
                          else result$predictor = eval(expr@formula[[3]], envir);
                  },			##   right hand formula
                  result=NULL); 
         },
         ###############################################
         ##  In function case, the user must enter the##
         ##  arguments of function into a list with   ##
         ##  names of function argument inside it, not##
         ##  not arguments seperatlly.                ##
         "function"={										##  function case                            ##
           .nf <- length(expr@fnc)       ##  length of all argments and function body.##
           .arglist <- expr@arguments		##  default arguments of the object in slot  ##
           ##                               arguments   ##
           temp <- formals(expr@fnc)   ##  list of argumetns and default values     ##
           .arglist[names(temp)] <- temp ##  A list of arguments and default values,  ##
           ## and object slot 'arguments' if not changed##
           .arglist[names(envir)] <- envir  ##  Substitude data to arguments, new values ##
           ##  will be replaced in default values and   ##
           ##  non default values authomatically replaced#
           ##  by called values.                        ##
           result <- do.call(expr@fnc,.arglist)
         }, 
         result=NULL
  )
  if(! is.atomic(result)){
    if(! is.null(result$response))                          ##  if any of the response or predictor be   ##
      if(any(is.nan(result$response))) result=Fault(FN=18)  ##  'nan' , NULL value return back.          ##
    if(! is.null(result$predictor))                         ##  due to any reason it can be, zero devide ##
      if(any(is.nan(result$predictor))) result=Fault(FN=18) ##  missing values ...., in this case feature##
  }
  return(result)															   ##  missing values ...., in this case feature##
}
setMethod("evald","nl.form",definition=eval.nl.form)
#############################################################|
##   Method $, nl.form.                                       \
###################################################            \
#method.skeleton("$","nl.form")                                 \
#getGeneric("$")                  # discover the generic method  |
#showMethods("$")                                                |
#---------------------------------------------------------------/
setMethod("$", "nl.form",
    function (x, name)  
	slot(x,name)
)
																			##  modification will be needed.             ##
#############################################################|
##   Method replacement, nl.form.                             \
###################################################            \
#--------------------------------------------------------------/
setReplaceMethod(f="[", signature="nl.form",
          definition=function (x,i,value){
            if (i=="selfStart"){x@selfStart<-value}
            validObject(x)
            return(x)
          }  
)
#same--------------------------------------------------------------/
setReplaceMethod(f="$", signature="nl.form",
                 definition=function (x,name,value){
                  if (name=="selfStart"){x@selfStart<-value}
                  else stop("There is no slot called",name)
                  return(x)
                 }  
)
# the above can be replaced by setMethod(f = "$<-", signature = "nl.form",
###################################################          |
##   Method selfStart, nl.form.                               \
##     ex: selfStart(nl.form object, initialize function)       \
##     (parameters is deleted here becuse it stored in object)   \
###################################################              |
setMethod("selfStart",
    signature(model = "nl.form"),
    function (model, initial, parameters, template) model@selfStart <- initial)

###################################################
##   Method getItial, nl.form.
##      ex:  getInitial(nl.form object, data,...)
###################################################
setMethod("getInitial","nl.form",
          function(object,data,...) {
            if(is.null(object@selfStart)) return(object$par)
            else return(object@selfStart(data,...) )
          }
          )

###################################################
##   Method all.vars, nl.form.
##      ex:  allvars(nl.form object,expr, 
##        functions = FALSE, max.names = -1L, unique = TRUE))
##        like original function, see help.
###################################################
#method.skeleton("all.vars","nl.form")
setMethod("all.vars",
          signature(expr = "nl.form"),
          function (expr, functions = FALSE, max.names = -1L, unique = TRUE) 
          {
            c(names(expr$par),expr$dependent,expr$independent,names(expr$arguments))
          }
)

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.form' &                                 |
#|                                                                                 |
#|                          eval  function                                         |
#|                                                                                 |
#|                            24/09/2008                                           |
#|                                                                                 |
#|                          selfStart & getInitial extended                        |
#|                                                                                 |
#|                               10/08/2011                                        |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, NSPEM                               |
#|                                                                                 |
#|                         Developed for R 02/10/2012                              |
#+#################################################################################+

					


