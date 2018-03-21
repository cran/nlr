#################################################################################|
#|        Object: robust measures. contains a masure, cut of points and plot      \
#|                                                                                |
#################################################################################/
setClass("nl.robmeas",representation(
											measure="numeric",
											cutofpoint = "numeric",
											name = "character"
											)
			)
#################################
##      constructor            ##
#################################
nl.robmeas<-function(measure,cutofpoint,name){
	new("nl.robmeas",measure=measure,cutofpoint=cutofpoint,name=name)
}

###################################################
##   Method $, nl.fitt.
###################################################
setMethod("$",
	signature(x = "nl.robmeas"),
	function(x,name){
		slot(x,name)
	} 
)

#######################################
##  plot method                      ##
##  (...) aruments to plot options   ##
#######################################
setMethod("plot",signature(x="nl.robmeas"),

#****************   plot a list of objects of any type   ****************
	function(x,y,...){
		dot.args <- list(...)
		if(is.null(dot.args$main)) main <- x@name
		else main <- dot.args$man
		if(is.null(dot.args$ylab)) ylab <- "Measure value"
		else ylab <- dot.args$ylab
		if(any(is.nan(x@measure))){
		  print("NAN in robust measure can not plot it these are the measure values")
      print(x@measure)
      print("cut of point of measure are equal:")
      print(x@cutofpoint)
		} 
    else {
      plot(x@measure,ylim=c(min(x@cutofpoint,x@measure),max(x@cutofpoint,x@measure)),ylab=ylab,main=main, 
           ...)
		  usr <- par("usr")
		  xlen <- usr[2] - usr[1]
		  x1 <- usr[1]+.03*xlen
		  x2 <- usr[2]-.03*xlen
		  for(i in x@cutofpoint)	lines(c(x1,x2),c(i,i))		
    }
	}
)

#+#################################################################################+
#|                                                                                 |
#|                   End of the object 'nl.robmeas'                                |
#|                                                                                 |
#|                             OCT 2008                                            |
#|                                                                                 |
#|                    Hossein Riazoshams, UPM, INSPEM                              |
#|                                                                                 |
#+#################################################################################+
