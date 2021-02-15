library(ggplot2)
library(dplyr)
library(gridExtra)
library(plotly)

depthGramPlot<-function(DG, limits=F, ids=NULL, print=F, plotly=F, plot.title="", col=NULL,shorten=T, pch=19, sp=2,st=4,sa=10, text.labels=""){  
  #Function to plot the 3 DepthGram representations from the output of the depthGram function
  
  ### DG: output of the depthGram function
  ### limits: Should empirical limits for outlier detection be drawn?
  ### ids: labels for data points
  ### print: Should the graphical output be optimized for printed version?
  ### plotly: Should the graphical output be displayed as an interactive plotly object?
  ### plot.title: main title for plot
  ### The remaining arguments are graphical parameters (text.labels is overridden by limits=T, for which only outliers labels are shown)
  
  n=length(DG$mei.mbd.d)
  
  type=c("Dimensions DepthGram","Time DepthGram","Time/Correlation DepthGram")
  
  if(is.null(ids)){
    ids=as.character(1:n)
  }
  
  DG=data.frame(ID=rep(ids,3),mei.mbd=c(DG$mei.mbd.d,DG$mei.mbd.t,DG$mei.mbd.t2), mbd.mei=c(DG$mbd.mei.d, DG$mbd.mei.t, DG$mbd.mei.t2), 
                type=rep(type,each=n))

  
  if(is.null(col)){
    hues = seq(15, 375, length=n+1)
    color <- hcl(h=hues, l=65, c=100)[1:n]
  }else{
    color=col
  }
  
  out=NULL
  if(limits){
    P2<-function(x,n){a0=2/n;a2=-n/(2*(n-1));return(a0+x+a2*x^2)}
    meis=seq(0,1,,n)
    DG$meis=rep(meis,3)
    DG$par=P2(DG$meis,n)
    distp=DG$mbd.mei-P2(1-DG$mei.mbd,n)
    q3=quantile(distp,0.75)
    q1=quantile(distp,0.25)
    DG$par2=DG$par+q3+1.5*(q3-q1)
    out<-unique(which(distp>q3+1.5*(q3-q1))%%n)
    pch=rep(1,n) #empty circle
    pch[out]=19 #solid circle 
    text.labels=rep("",n) 
    if (!plotly){
      if (shorten){ #shorten text labels
        text.labels[out]=sapply(ids[out],function(x) substr(x,1,min(15,nchar(x))))
      }else{
        text.labels[out]=ids[out]
      }
    }
    if(is.null(col)){
      hues = seq(15, 375, length=length(out)+1)
      color.out <- hcl(h=hues, l=65, c=100)[1:length(out)]
      color<-rep(8,n)
      color[out]<-color.out
    }
  }
  
  if(print){
    sp=3
    st=5
    sa=14
  }
  
  plots=list()
  
  for (i in 1:3){
    
    dat=DG[which(DG$type==type[i]),]
    
    plots[[i]]<- ggplot(dat, aes(x=1-mei.mbd, y=mbd.mei)) + 
      geom_point(aes(x=1-mei.mbd, y=mbd.mei, group=ID),color=color,size=sp,shape=pch) +
      geom_text(aes(x=1-mei.mbd, y=mbd.mei),label=text.labels,color=color,hjust=-0.15, vjust=-0.15,size=st)+
      xlim(c(0,1.005))+ ylim(c(0,0.525)) +theme_minimal()+
      theme(axis.title= element_text(size=sa),axis.text= element_text(size=sa-2),title=element_text(size=sa))
    
    if(limits==TRUE){
      plots[[i]]<- plots[[i]] +geom_line(aes(x=meis, y=par),col=1,na.rm=T)+
        geom_line(aes(x=meis, y=par2),col=1,lty=2,na.rm=T)
    }
    
    if (i==1){pt<-plot.title}else{pt<-""}
    
    if (plotly){
      plots[[i]]<-plots[[i]]+facet_wrap(~type)+ggtitle(pt) 
    }else{
      if(i<=2){
        plots[[i]]<-plots[[i]]+labs(title=pt,subtitle=type[i],x=expression("1-MEI("~paste(bold('MBD'[d])~")")), 
                                    y=expression("MBD("~paste(bold('MEI'[d])~")")))
      }else{
        plots[[i]]<-plots[[i]] + labs(title=pt,subtitle=type[i], x=expression("1-MEI("~paste(bold(widetilde(MBD)[t])~")")),
                                      y=expression("MBD("~paste(bold(widetilde(MEI)[t])~")")))
      }
    }
    
  }
  
  if(plotly){
    p<-subplot(plots,shareY = TRUE,shareX =TRUE)%>%
      layout(title = plot.title, yaxis = list(title="MBD(MEI)"), xaxis = list(title="1-MEI(MBD)"))
    print(p)
  }else{
    p<-do.call("grid.arrange",c(plots,nrow=1))
  }

  return(list(p=list(dimDG = plots[[1]], timeDG = plots[[2]], corrDG = plots[[3]], fullDG = p), out=out, color=color))
}
