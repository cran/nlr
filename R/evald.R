evald <-function (expr, envir = parent.frame(), enclos = if (is.list(envir) || 
    is.pairlist(envir)) parent.frame() else baseenv(),...) {
	lst<-list(...)
	if(! is.environment(envir))
		if(length(lst)!=0){
			
			lst2<-c(envir,lst)       # it will be a list, even if envir is a data.frame
			return(eval(expr,envir=lst2))
		}
		else return(eval(expr,envir))
	else
		if(length(lst)!=0){
			lst2<-c(envir,lst)
			return(eval(expr))#,envir=lst2)
		}
		else 
			return(eval(expr))
}

#a=y~x
#eval(a,list(x=3,y=4))
#eval(a)
#b2=evald(a)
#b2
#eval(a[[3]])
#evald(a[[3]],b=4)

#b=function(a,...)
#evald(a,...)
#b2=b(a[[3]],list(x=4),y=5)
#b2
#b2=evald(a,y=5)
#b2

