#+#########################################################################+
###    method.db data base of fitted methods                               |
#+#########################################################################+

#+#########################################################################+
###    Method object. store computed methods                                \
###       Slots:                                                            /
###             Method in program, depends on methods contains:            /
###                   Method: character of method name                     \
###                   Detail1: detail of the method name                    \ 
###                   Detail2: deatil of computation method                 |
###                   MethodID: ID of the method                           /
###                   Subroutine: the function when the calculation is done\
###                                                                         |
###                                                                         |
############################################################################+


setClass("fittmethod",representation(
										methodID=           "numeric",
										method=             "character",
										detail=             "character",
										methodBR=           "numeric",           ## method branch like computation method
										detailBR=           "character",
										subroutine=         "character",
										lossfunction =      "character",
										subroutineBR =      "character"
										)
			)
setClassUnion("fittmethodorNULL", c("fittmethod", "NULL"))
#################################
##      constructor            ##
#################################
fittmethod <- function(methodID,
						method=		(nlr::db.method$method[nlr::db.method$methodID==methodID]),
						detail=		as.character(nlr::db.method$detail[nlr::db.method$methodID==methodID]),
						methodBR= 		(nlr::db.method$methodBR[nlr::db.method$methodID==methodID]),
						detailBR=		as.character(nlr::db.methodBR$detailBR[nlr::db.methodBR$methodBR==methodBR]),
						subroutine=	as.character(nlr::db.method$subroutine[nlr::db.method$methodID==methodID]),
						lossfunction= "",
						subroutineBR = as.character(nlr::db.methodBR$subroutineBR[nlr::db.methodBR$methodBR==methodBR])
					){
	return( new( "fittmethod",methodID=methodID,method=method,detail=detail, methodBR=methodBR,detailBR=detailBR,
		subroutine=subroutine,lossfunction=lossfunction,subroutineBR=subroutineBR))
}

###################################################
##   fittmethod $, Fault.
###################################################
setMethod("$","fittmethod",
	function(x,name){
		slot(x,name)
	} 
)
###################################################
##   fittmethod summary.
setMethod("summary",
          signature(object = "fittmethod"),
          function (object, ...) 
          {
            cat("-------------------------------------------------------------------------------------------", "\n")
            cat("methodID , method , detail , methodBR , detailBR , subroutine , lossfunction , subroutineBR", "\n")
            cat("-------------------------------------------------------------------------------------------", "\n")
            cat(object$methodID,object$method,object$detail,object$methodBR,object$detailBR,object$subroutine,
                object$lossfunction,object$subroutineBR,"\n")
            cat("-------------------------------------------------------------------------------------------", "\n")
          }
)

###################################################
##   fittmethod print.
setMethod("print",
          signature(x = "fittmethod"),
          function (x, ...) 
          {
            cat("-------------------------------------------------------------------------------------------", "\n")
            cat("methodID , method , detail , methodBR , detailBR , subroutine , lossfunction , subroutineBR", "\n")
            cat("-------------------------------------------------------------------------------------------", "\n")
            cat(x$methodID,x$method,x$detail,x$methodBR,x$detailBR,x$subroutine,
                x$lossfunction,x$subroutineBR,"\n")
            cat("-------------------------------------------------------------------------------------------", "\n")
            
          }
)

#+#################################################################+
#|                   End of the object 'method'                    |
#|                            03/12/2012                           |
#|                    Hossein Riazoshams, Statistics Department    |
#|                      stockholm University                       |
#+#################################################################+
# fittmethod(1,"ss","ss",1,"ss","ss")
#fittmethod(1)
#removeClass("fittmethod")
#removeMethod("$","fittmethod")
#remove("fittmethod")