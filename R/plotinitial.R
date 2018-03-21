plotinitial<-function(form,data,start=getInitial(form,data),length.out=100,...){
  datalist <-c(start,data)
  xd <- unlist(data[form$independent])
  n<-length(xd)
  xsc<-seq(min(xd),max(xd),length.out=length.out)
  datalist[[form$independent]]=xsc
  ysc<-as.numeric(eval(form,datalist)$predictor)

  x1<-min(xsc,data[[form$independent]])
  x2<-max(xsc,data[[form$independent]])
  y1<-min(ysc,data[[form$dependent]])
  y2<-max(ysc,data[[form$dependent]])
  x<-data[[form$independent]]
  y<-data[[form$dependent]]
  plot(x,y,
              xlim=c(x1,x2),
              ylim=c(y1,y2),
              main=paste("Initial values \n",form@name,"model")
       )
  lines(xsc,ysc)
}


plotinitial2<-function(form,data,start,length.out=100){
  datalist <-c(start,data)
  xd <- unlist(data[form$independent])
  n<-length(xd)
  xsc<-seq(min(xd),max(xd),length.out=length.out)
  datalist[[form$independent]]=xsc
  ysc<-as.numeric(eval(form,datalist)$predictor)
  
  x1<-min(xsc,data[[form$independent]])
  x2<-max(xsc,data[[form$independent]])
  y1<-min(ysc,data[[form$dependent]])
  y2<-max(ysc,data[[form$dependent]])
  x<-data[[form$independent]]
  y<-data[[form$dependent]]
  plot(xsc,ysc,type="l")
  plot(x,y,
       xlim=c(x1,x2),
       ylim=c(y1,y2))
  lines(xsc,ysc)
}


remove(plotinitial2)
